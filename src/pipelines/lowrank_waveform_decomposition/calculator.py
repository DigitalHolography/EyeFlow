"""Pure NumPy low-rank decomposition of beat-aligned segment waveforms."""

from __future__ import annotations

import numpy as np


class LowRankWaveformDecompositionCalculator:
    """Compute joint and beat-local SVD endpoint panels.

    Input waveforms use EyeFlow's native ``(sample, beat, branch, radius)``
    order.  The joint decomposition treats each beat/location waveform as one
    column of a ``(sample, beat * branch * radius)`` matrix after subtracting
    the temporal mean of every column.
    """

    eps = 1e-12
    min_valid_samples_fraction = 0.95
    min_valid_columns = 3
    exported_modes = 2

    def compute(self, waveforms: np.ndarray, beat_period_seconds: np.ndarray) -> dict:
        periods = self._normalize_periods(beat_period_seconds)
        block = self._ensure_segment_shape(waveforms, periods.size)
        return self._compute_joint(block, periods)

    def compute_per_beat(
        self,
        waveforms: np.ndarray,
        beat_period_seconds: np.ndarray,
    ) -> dict[str, np.ndarray]:
        """Run an independent ``(sample, branch * radius)`` SVD per beat."""
        periods = self._normalize_periods(beat_period_seconds)
        block = self._ensure_segment_shape(waveforms, periods.size)
        n_beats = block.shape[1]
        names = (
            "A1_b_pb",
            "A2_b_pb",
            "R1_b_pb",
            "R2_b_pb",
            "rho1_b_pb",
            "rho2_b_pb",
            "R0_b_pb",
            "TPR_b_pb",
            "rho0_b_pb",
            "MPR_b_pb",
            "effective_rank_b_pb",
            "participation_ratio_b_pb",
            "alpha_b_pb",
            "G1_b_pb",
        )
        result = {
            name: np.full((n_beats,), np.nan, dtype=float) for name in names
        }

        for beat in range(n_beats):
            if not np.isfinite(periods[beat]) or periods[beat] <= 0:
                continue
            panel = self._compute_single_beat(block[:, beat], periods[beat])
            for name in names:
                result[name][beat] = panel[name.removesuffix("_b_pb")]
        return result

    @staticmethod
    def _normalize_periods(periods: np.ndarray) -> np.ndarray:
        values = np.asarray(periods, dtype=float)
        if values.ndim == 1:
            return values
        if values.ndim == 2 and 1 in values.shape:
            return values.reshape(-1)
        raise ValueError(
            "Beat periods must have shape (beat,), (1, beat), or (beat, 1); "
            f"got {values.shape}."
        )

    @staticmethod
    def _ensure_segment_shape(block: np.ndarray, n_beats: int) -> np.ndarray:
        values = np.asarray(block, dtype=float)
        if values.ndim != 4:
            raise ValueError(
                "Segment waveforms must have shape "
                "(sample, beat, branch, radius); "
                f"got {values.shape}."
            )
        if values.shape[1] == n_beats:
            return values
        if values.shape[0] == n_beats and values.shape[1] != n_beats:
            return np.transpose(values, (1, 0, 2, 3))
        raise ValueError(
            "No waveform axis matches the beat-period count "
            f"({n_beats}) in shape {values.shape}."
        )

    def _compute_joint(self, block: np.ndarray, periods: np.ndarray) -> dict:
        n_t, n_beats, n_branches, n_radii = block.shape
        finite_fraction = np.mean(np.isfinite(block), axis=0)
        period_valid = np.isfinite(periods) & (periods > 0)
        valid = finite_fraction >= self.min_valid_samples_fraction
        valid &= period_valid[:, None, None]

        mu = self._column_mean(block)
        centered = np.where(np.isfinite(block), block, mu[None]) - mu[None]
        centered = np.where(np.isfinite(centered), centered, 0.0)
        total_rms = np.sqrt(np.mean(centered**2, axis=0))
        total_rms[~valid] = np.nan

        valid_counts = np.sum(valid, axis=(1, 2))
        out: dict[str, object] = {
            "shape": {
                "n_t": n_t,
                "n_beats": n_beats,
                "n_branches": n_branches,
                "n_radii": n_radii,
                "n_total_columns": int(n_beats * n_branches * n_radii),
                "n_valid_columns": int(np.sum(valid)),
            },
            "finite_fraction_per_column": finite_fraction,
            "valid_column_mask": valid,
            "valid_counts_per_beat": valid_counts,
            "valid_fraction_per_beat": valid_counts
            / float(max(1, n_branches * n_radii)),
            "beat_period_valid": period_valid,
            "mu": mu,
            "x_full": centered,
            "total_rms_bkr": total_rms,
        }
        local_mpr = self._safe_ratio(np.abs(mu), total_rms)
        local_mpr[~valid] = np.nan
        beatwise, acq = self._baseline_endpoints(
            mu,
            total_rms,
            local_mpr,
            valid,
            periods,
        )
        out["mean_to_pulsatile_ratio_bkr"] = local_mpr
        out["mean_pulsatile_ratio_bkr"] = local_mpr
        out["beatwise"] = beatwise
        out["acq"] = acq

        if int(np.sum(valid)) < self.min_valid_columns:
            out["svd_available"] = False
            out["svd_reason"] = "too_few_valid_columns"
            return out

        matrix = centered[:, valid]
        if matrix.size == 0:
            out["svd_available"] = False
            out["svd_reason"] = "empty_valid_matrix"
            return out

        U, singular_values, Vt = np.linalg.svd(matrix, full_matrices=False)
        energy = singular_values**2
        total_energy = float(np.sum(energy))
        energy_fraction = (
            energy / total_energy
            if total_energy > self.eps
            else np.full(energy.shape, np.nan, dtype=float)
        )
        n_modes = min(self.exported_modes, singular_values.size)
        scores: list[np.ndarray] = []
        sign_flips = np.zeros((self.exported_modes,), dtype=np.int32)

        for mode in range(n_modes):
            mode_scores = singular_values[mode] * Vt[mode]
            if self._safe_nanmedian(mode_scores) < 0:
                U[:, mode] *= -1.0
                Vt[mode] *= -1.0
                mode_scores *= -1.0
                sign_flips[mode] = 1
            scores.append(mode_scores)

        u_panel = np.full((n_t, self.exported_modes), np.nan, dtype=float)
        score_panel = np.full(
            (self.exported_modes, n_beats, n_branches, n_radii),
            np.nan,
            dtype=float,
        )
        rms_mode_panel = np.full_like(score_panel, np.nan)
        residual_rms_panel = np.full_like(score_panel, np.nan)
        residual_t_panel = np.full(
            (self.exported_modes, n_t, n_beats, n_branches, n_radii),
            np.nan,
            dtype=float,
        )

        for mode in range(n_modes):
            u_panel[:, mode] = U[:, mode]
            score_panel[mode, valid] = scores[mode]
            rms_u = float(np.sqrt(np.mean(U[:, mode] ** 2)))
            rms_mode_panel[mode, valid] = np.abs(scores[mode]) * rms_u

            reconstruction = U[:, : mode + 1] @ np.vstack(scores[: mode + 1])
            residual = matrix - reconstruction
            residual_rms_panel[mode, valid] = np.sqrt(np.mean(residual**2, axis=0))
            residual_t_panel[mode][:, valid] = residual

            index = mode + 1
            self._append_mode_endpoints(
                beatwise,
                acq,
                index,
                rms_mode_panel[mode],
                residual_rms_panel[mode],
                score_panel[mode],
                total_rms,
                valid,
            )

        acq.update(self._spectrum_endpoints(singular_values, energy_fraction))
        out.update(
            {
                "svd_available": True,
                "svd_reason": "ok",
                "n_modes_panel": n_modes,
                "U_panel": u_panel,
                "score_panel_bkr": score_panel,
                "rms_mode_panel": rms_mode_panel,
                "residual_rms_panel": residual_rms_panel,
                "residual_t_bkr_panel": residual_t_panel,
                "sign_flips": sign_flips,
                "s": singular_values,
                "energy": energy,
                "energy_fraction": energy_fraction,
            }
        )
        return out

    def _baseline_endpoints(
        self,
        mu: np.ndarray,
        total_rms: np.ndarray,
        local_mpr: np.ndarray,
        valid: np.ndarray,
        periods: np.ndarray,
    ) -> tuple[dict, dict]:
        mu_b = self._median_per_beat(mu, valid)
        abs_mu_b = self._median_per_beat(np.abs(mu), valid)
        r0_b = self._median_per_beat(total_rms, valid)
        mpr_b = self._median_per_beat(local_mpr, valid)
        rho0_b = np.where(np.isfinite(r0_b) & (r0_b > self.eps), 1.0, np.nan)
        r0 = self._safe_nanmedian(total_rms[valid])
        beatwise = {
            "mu_b": mu_b,
            "abs_mu_b": abs_mu_b,
            "R0_b": r0_b,
            # Compatibility alias retained for the first AngioEye wrapper.
            "TPR_b": r0_b,
            "rho0_b": rho0_b,
            "MPR_b": mpr_b,
            "mpr_b": mpr_b,
            "spatial_mad_mu_b": self._spatial_mad_per_beat(mu, valid),
        }
        acq = {
            "mu_acq": self._safe_nanmedian(mu[valid]),
            "abs_mu_acq": self._safe_nanmedian(np.abs(mu[valid])),
            "beat_period_mean": self._safe_nanmean(periods[periods > 0]),
            "beat_period_median": self._safe_nanmedian(periods[periods > 0]),
            "beat_period_std": self._safe_nanstd(periods[periods > 0]),
            "R0": r0,
            # Compatibility alias retained for the first AngioEye wrapper.
            "TPR": r0,
            "rho0": self._safe_nanmedian(rho0_b),
            "MPR": self._safe_nanmedian(local_mpr[valid]),
            "sigma_mu_beat": self._safe_nanstd(mu_b),
            "mad_mu_beat": self._safe_nanmad(mu_b),
            "sigma_R0_beat": self._safe_nanstd(r0_b),
            "mad_R0_beat": self._safe_nanmad(r0_b),
            "sigma_rho0_beat": self._safe_nanstd(rho0_b),
            "mad_rho0_beat": self._safe_nanmad(rho0_b),
            "cv_rho0_beat": self._safe_nancv(rho0_b),
            "spatial_mad_mu_median_over_beats": self._safe_nanmedian(
                beatwise["spatial_mad_mu_b"]
            ),
        }
        acq["mpr"] = acq["MPR"]
        acq["sigma_TPR_beat"] = acq["sigma_R0_beat"]
        acq["mad_TPR_beat"] = acq["mad_R0_beat"]
        acq["sigma_mpr_beat"] = self._safe_nanstd(mpr_b)
        acq["mad_mpr_beat"] = self._safe_nanmad(mpr_b)
        acq["cv_mpr_beat"] = self._safe_nancv(mpr_b)
        return beatwise, acq

    def _append_mode_endpoints(
        self,
        beatwise: dict,
        acq: dict,
        mode: int,
        amplitude: np.ndarray,
        residual: np.ndarray,
        scores: np.ndarray,
        total_rms: np.ndarray,
        valid: np.ndarray,
    ) -> None:
        amplitude_b = self._median_per_beat(amplitude, valid)
        residual_b = self._median_per_beat(residual, valid)
        rho_b = self._safe_ratio(residual_b, beatwise["R0_b"])
        median_score_b = self._median_per_beat(np.abs(scores), valid)
        spatial_mad_a_b = self._spatial_mad_per_beat(amplitude, valid)
        spatial_mad_r_b = self._spatial_mad_per_beat(residual, valid)

        beatwise[f"A{mode}_b"] = amplitude_b
        beatwise[f"R{mode}_b"] = residual_b
        beatwise[f"rho{mode}_b"] = rho_b
        beatwise[f"median_abs_a{mode}_b"] = median_score_b
        beatwise[f"spatial_mad_A{mode}_b"] = spatial_mad_a_b
        beatwise[f"spatial_mad_R{mode}_b"] = spatial_mad_r_b

        acq[f"A{mode}"] = self._safe_nanmedian(amplitude[valid])
        acq[f"R{mode}"] = self._safe_nanmedian(residual[valid])
        acq[f"rho{mode}"] = self._safe_scalar_ratio(
            acq[f"R{mode}"], acq["R0"]
        )
        acq[f"median_abs_a{mode}"] = self._safe_nanmedian(np.abs(scores)[valid])
        for prefix, values in (
            (f"A{mode}", amplitude_b),
            (f"R{mode}", residual_b),
            (f"rho{mode}", rho_b),
        ):
            acq[f"sigma_{prefix}_beat"] = self._safe_nanstd(values)
            acq[f"mad_{prefix}_beat"] = self._safe_nanmad(values)
            acq[f"cv_{prefix}_beat"] = self._safe_nancv(values)
        acq[f"spatial_mad_A{mode}_median_over_beats"] = self._safe_nanmedian(
            spatial_mad_a_b
        )
        acq[f"spatial_mad_R{mode}_median_over_beats"] = self._safe_nanmedian(
            spatial_mad_r_b
        )

    def _compute_single_beat(self, beat: np.ndarray, period: float) -> dict[str, float]:
        finite_fraction = np.mean(np.isfinite(beat), axis=0)
        valid = finite_fraction >= self.min_valid_samples_fraction
        mu = self._column_mean(beat)
        centered = np.where(np.isfinite(beat), beat, mu[None]) - mu[None]
        centered = np.where(np.isfinite(centered), centered, 0.0)
        total_rms = np.sqrt(np.mean(centered**2, axis=0))
        r0 = self._safe_nanmedian(total_rms[valid])
        local_mpr = self._safe_ratio(np.abs(mu), total_rms)
        empty = {
            "A1": np.nan,
            "A2": np.nan,
            "R1": np.nan,
            "R2": np.nan,
            "rho1": np.nan,
            "rho2": np.nan,
            "R0": r0,
            "TPR": r0,
            "rho0": 1.0 if np.isfinite(r0) and r0 > self.eps else np.nan,
            "MPR": self._safe_nanmedian(local_mpr[valid]),
            "effective_rank": np.nan,
            "participation_ratio": np.nan,
            "alpha": np.nan,
            "G1": np.nan,
        }
        if int(np.sum(valid)) < self.min_valid_columns:
            return empty

        matrix = centered[:, valid]
        U, singular_values, Vt = np.linalg.svd(matrix, full_matrices=False)
        energy = singular_values**2
        total_energy = float(np.sum(energy))
        energy_fraction = (
            energy / total_energy
            if total_energy > self.eps
            else np.full(energy.shape, np.nan, dtype=float)
        )
        empty.update(self._spectrum_endpoints(singular_values, energy_fraction))
        empty.pop("eta1", None)
        empty.pop("eta2", None)
        empty.pop("eta12", None)

        scores: list[np.ndarray] = []
        for mode in range(min(self.exported_modes, singular_values.size)):
            mode_scores = singular_values[mode] * Vt[mode]
            scores.append(mode_scores)
            rms_u = float(np.sqrt(np.mean(U[:, mode] ** 2)))
            amplitude = np.abs(mode_scores) * rms_u
            residual = matrix - U[:, : mode + 1] @ np.vstack(scores)
            residual_rms = np.sqrt(np.mean(residual**2, axis=0))
            index = mode + 1
            empty[f"A{index}"] = self._safe_nanmedian(amplitude)
            empty[f"R{index}"] = self._safe_nanmedian(residual_rms)
            empty[f"rho{index}"] = self._safe_scalar_ratio(
                empty[f"R{index}"], r0
            )
        return empty

    def _spectrum_endpoints(
        self,
        singular_values: np.ndarray,
        energy_fraction: np.ndarray,
    ) -> dict[str, float]:
        energy = np.asarray(singular_values, dtype=float) ** 2
        eta1 = float(energy_fraction[0]) if energy_fraction.size >= 1 else np.nan
        eta2 = float(energy_fraction[1]) if energy_fraction.size >= 2 else np.nan
        leading_energy = float(energy[0]) if energy.size >= 1 else np.nan
        alpha = self._safe_scalar_ratio(
            float(np.sum(energy[1:])), leading_energy
        )
        g1 = (
            self._safe_scalar_ratio(
                float(singular_values[0] - singular_values[1]),
                float(singular_values[0]),
            )
            if singular_values.size >= 2
            else np.nan
        )
        return {
            "eta1": eta1,
            "eta2": eta2,
            "eta12": float(np.sum(energy_fraction[:2])),
            "spectrum_mode_count": int(energy_fraction.size),
            "effective_rank": self._effective_rank(energy_fraction),
            "participation_ratio": self._participation_ratio(energy_fraction),
            "alpha": alpha,
            "G1": g1,
        }

    @staticmethod
    def _column_mean(values: np.ndarray) -> np.ndarray:
        finite = np.isfinite(values)
        counts = np.sum(finite, axis=0)
        sums = np.sum(np.where(finite, values, 0.0), axis=0)
        return np.divide(
            sums,
            counts,
            out=np.full_like(sums, np.nan, dtype=float),
            where=counts > 0,
        )

    def _median_per_beat(self, values: np.ndarray, valid: np.ndarray) -> np.ndarray:
        result = np.full((values.shape[0],), np.nan, dtype=float)
        for beat in range(values.shape[0]):
            result[beat] = self._safe_nanmedian(values[beat][valid[beat]])
        return result

    def _spatial_mad_per_beat(
        self,
        values: np.ndarray,
        valid: np.ndarray,
    ) -> np.ndarray:
        result = np.full((values.shape[0],), np.nan, dtype=float)
        for beat in range(values.shape[0]):
            result[beat] = self._safe_nanmad(values[beat][valid[beat]])
        return result

    def _effective_rank(self, fractions: np.ndarray) -> float:
        p = np.asarray(fractions, dtype=float)
        p = p[np.isfinite(p) & (p > 0)]
        return float(np.exp(-np.sum(p * np.log(p)))) if p.size else np.nan

    @staticmethod
    def _participation_ratio(fractions: np.ndarray) -> float:
        p = np.asarray(fractions, dtype=float)
        p = p[np.isfinite(p) & (p > 0)]
        denominator = float(np.sum(p**2))
        return float(1.0 / denominator) if denominator > 0 else np.nan

    def _safe_ratio(self, numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
        num, den = np.broadcast_arrays(
            np.asarray(numerator, dtype=float), np.asarray(denominator, dtype=float)
        )
        return np.divide(
            num,
            den,
            out=np.full(num.shape, np.nan, dtype=float),
            where=np.isfinite(num) & np.isfinite(den) & (np.abs(den) > self.eps),
        )

    def _safe_scalar_ratio(self, numerator: float, denominator: float) -> float:
        if (
            not np.isfinite(numerator)
            or not np.isfinite(denominator)
            or abs(denominator) <= self.eps
        ):
            return np.nan
        return float(numerator / denominator)

    @staticmethod
    def _safe_nanmean(values: np.ndarray) -> float:
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        return float(np.mean(finite)) if finite.size else np.nan

    @staticmethod
    def _safe_nanmedian(values: np.ndarray) -> float:
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        return float(np.median(finite)) if finite.size else np.nan

    @staticmethod
    def _safe_nanstd(values: np.ndarray) -> float:
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        return float(np.std(finite)) if finite.size else np.nan

    @classmethod
    def _safe_nanmad(cls, values: np.ndarray) -> float:
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        if not finite.size:
            return np.nan
        median = np.median(finite)
        return float(np.median(np.abs(finite - median)))

    def _safe_nancv(self, values: np.ndarray) -> float:
        mean = self._safe_nanmean(values)
        std = self._safe_nanstd(values)
        return self._safe_scalar_ratio(std, abs(mean))


__all__ = ["LowRankWaveformDecompositionCalculator"]
