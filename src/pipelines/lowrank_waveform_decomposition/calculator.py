"""Low-rank decomposition of beat-aligned segment waveforms."""

from __future__ import annotations

import warnings

import numpy as np


def normalize_periods(periods: np.ndarray) -> np.ndarray:
    """Return beat periods as shape (1, n_beats)."""
    values = np.asarray(periods, dtype=float)
    if values.ndim == 1:
        return values.reshape(1, -1)
    if values.ndim == 2 and values.shape[0] == 1:
        return values
    if values.ndim == 2 and values.shape[1] == 1:
        return values.T
    raise ValueError(
        "Beat periods must have shape (beat,), (1, beat), or (beat, 1); "
        f"got {values.shape}."
    )


def ensure_segment_shape(
    block: np.ndarray, periods: np.ndarray | None = None
) -> np.ndarray:
    """Ensure waveforms are (sample, beat, branch, radius)."""
    values = np.asarray(block, dtype=float)
    if values.ndim != 4:
        raise ValueError(
            "Segment waveforms must have shape "
            f"(sample, beat, branch, radius); got {values.shape}."
        )
    if periods is None:
        return values
    n_beats = int(normalize_periods(periods).shape[1])
    if values.shape[1] == n_beats:
        return values
    if values.shape[0] == n_beats and values.shape[1] != n_beats:
        return np.transpose(values, (1, 0, 2, 3))
    raise ValueError(
        "No waveform axis matches the beat-period count "
        f"({n_beats}) in shape {values.shape}."
    )


def mean_subtract(v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Column-wise temporal mean subtraction for SVD input."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="Mean of empty slice", category=RuntimeWarning
        )
        mu = np.nanmean(v, axis=0)
    filled = np.where(np.isfinite(v), v, mu[None])
    centered = filled - mu[None]
    centered = np.where(np.isfinite(centered), centered, 0.0)
    return mu, centered


class LowRankWaveformDecompositionCalculator:
    """Compute joint and beat-local SVD endpoint panels.

    Input waveforms use EyeFlow\'s native ``(sample, beat, branch, radius)``
    order. Formulas follow AngioEye\'s manuscript-aligned definitions
    (two-stage beat medians for acquisition TPR/A/rho; alpha = (1-p1)/p1;
    G1 = 1 - lambda2/lambda1; mpr_prime = abs_mu_acq / TPR).
    """

    eps = 1e-12
    min_valid_samples_fraction = 0.95
    min_valid_columns = 3
    exported_modes = 2

    def compute(self, waveforms: np.ndarray, beat_period_seconds: np.ndarray) -> dict:
        """Joint SVD representation dict (EyeFlow / AngioEye shared shape)."""
        periods = normalize_periods(beat_period_seconds)
        block = ensure_segment_shape(waveforms, periods)
        return self._compute_representation(block, periods)

    def compute_per_beat(
        self,
        waveforms: np.ndarray,
        beat_period_seconds: np.ndarray,
    ) -> dict[str, np.ndarray]:
        """Independent per-beat SVD endpoint sequences."""
        periods = normalize_periods(beat_period_seconds)
        block = ensure_segment_shape(waveforms, periods)
        panels = self.per_beat_svd_panels(block, periods)
        return self._compute_per_beat_endpoints(panels)

    def compute_acquisition(
        self,
        waveforms: np.ndarray,
        beat_period_seconds: np.ndarray,
    ) -> dict | None:
        """Convenience bundle used by AngioEye\'s cohort adapter."""
        periods = normalize_periods(beat_period_seconds)
        block = ensure_segment_shape(waveforms, periods)
        rep = self._compute_representation(block, periods)
        if not rep.get("svd_available", False):
            # Still return baseline-only payload when SVD is unavailable so
            # callers can decide; mirror AngioEye\'s None only when no usable
            # waveform exists at all.
            if "acq" not in rep:
                return None
        per_beat = self._compute_per_beat_endpoints(
            self.per_beat_svd_panels(block, periods)
        )
        beat_period_arr = np.asarray(periods[0], dtype=float)
        return {
            "acq": rep["acq"],
            "beatwise": rep["beatwise"],
            "per_beat_svd": per_beat,
            "mu": rep["mu"],
            "energy_fraction": np.asarray(
                rep.get("energy_fraction", []), dtype=float
            ),
            "beat_period_mean": (
                float(np.nanmean(beat_period_arr))
                if beat_period_arr.size
                else float("nan")
            ),
            "beat_period_sd": (
                float(np.nanstd(beat_period_arr, ddof=1))
                if beat_period_arr.size > 1
                else float("nan")
            ),
            "beat_period_b": beat_period_arr,
            "valid_fraction_per_beat": np.asarray(
                rep.get("valid_fraction_per_beat", []), dtype=float
            ),
            "n_valid_columns": int(rep["shape"]["n_valid_columns"]),
            "n_total_columns": int(rep["shape"]["n_total_columns"]),
            "representation": rep,
        }

    # ------------------------------------------------------------------
    # Safe reducers / aggregation
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_nanmean(x: np.ndarray) -> float:
        x = np.asarray(x, dtype=float)
        if x.size == 0 or not np.any(np.isfinite(x)):
            return np.nan
        return float(np.nanmean(x))

    @staticmethod
    def _safe_nanmedian(x: np.ndarray) -> float:
        x = np.asarray(x, dtype=float)
        if x.size == 0 or not np.any(np.isfinite(x)):
            return np.nan
        return float(np.nanmedian(x))

    @staticmethod
    def _safe_nanstd(x: np.ndarray) -> float:
        x = np.asarray(x, dtype=float)
        if x.size == 0 or not np.any(np.isfinite(x)):
            return np.nan
        return float(np.nanstd(x))

    @staticmethod
    def _safe_nanmad(x: np.ndarray) -> float:
        x = np.asarray(x, dtype=float)
        if x.size == 0 or not np.any(np.isfinite(x)):
            return np.nan
        med = np.nanmedian(x)
        return float(np.nanmedian(np.abs(x - med)))

    def _safe_nancv(self, x: np.ndarray) -> float:
        mu = self._safe_nanmean(x)
        sd = self._safe_nanstd(x)
        if (not np.isfinite(mu)) or (not np.isfinite(sd)) or abs(mu) <= self.eps:
            return np.nan
        return float(sd / (abs(mu) + self.eps))

    def _median_kr_per_beat(
        self, arr_bkr: np.ndarray, valid_mask: np.ndarray
    ) -> np.ndarray:
        n_beats = int(arr_bkr.shape[0])
        out = np.full((n_beats,), np.nan, dtype=float)
        for b in range(n_beats):
            vals = np.asarray(arr_bkr[b], dtype=float)
            mask = np.asarray(valid_mask[b], dtype=bool)
            if not np.any(mask):
                continue
            x = vals[mask]
            if x.size == 0 or not np.any(np.isfinite(x)):
                continue
            out[b] = float(np.nanmedian(x))
        return out

    def _spatial_mad_per_beat(
        self, arr_bkr: np.ndarray, valid_mask: np.ndarray
    ) -> np.ndarray:
        n_beats = int(arr_bkr.shape[0])
        out = np.full((n_beats,), np.nan, dtype=float)
        for b in range(n_beats):
            vals = np.asarray(arr_bkr[b], dtype=float)
            mask = np.asarray(valid_mask[b], dtype=bool)
            if not np.any(mask):
                continue
            x = vals[mask]
            if x.size == 0 or not np.any(np.isfinite(x)):
                continue
            med = np.nanmedian(x)
            out[b] = float(np.nanmedian(np.abs(x - med)))
        return out

    def aggregate_beatwise(self, values_per_beat: np.ndarray, stat: str) -> float:
        x = np.asarray(values_per_beat, dtype=float)
        x = x[np.isfinite(x)]
        if x.size == 0:
            return float("nan")
        return float(np.mean(x)) if stat == "mean" else float(np.median(x))

    def aggregate_rho(self, R_b: np.ndarray, TPR_b: np.ndarray, stat: str) -> float:
        R = self.aggregate_beatwise(R_b, stat)
        T = self.aggregate_beatwise(TPR_b, stat)
        if not np.isfinite(R) or not np.isfinite(T) or T <= self.eps:
            return float("nan")
        return float(R / (T + self.eps))

    def _effective_rank(self, energy_fraction: np.ndarray) -> float:
        p = np.asarray(energy_fraction, dtype=float)
        p = p[np.isfinite(p) & (p > 0)]
        if p.size == 0:
            return np.nan
        return float(np.exp(-np.sum(p * np.log(p + self.eps))))

    def _participation_ratio(self, energy_fraction: np.ndarray) -> float:
        p = np.asarray(energy_fraction, dtype=float)
        p = p[np.isfinite(p) & (p > 0)]
        if p.size == 0:
            return np.nan
        denom = float(np.sum(p**2))
        if denom <= 0:
            return np.nan
        return float(1.0 / denom)

    def _alpha(self, s: np.ndarray, n_modes: int) -> float:
        """AngioEye formula: (1 - p1) / p1 on the leading exported modes."""
        s = np.asarray(s, dtype=float)
        if n_modes < 2 or s.size < 2:
            return np.nan
        energy_panel = s[:n_modes] ** 2
        p1 = float(energy_panel[0] / (np.sum(energy_panel) + self.eps))
        if not np.isfinite(p1) or p1 <= 0:
            return np.nan
        return float((1.0 - p1) / p1)

    def _g1(self, s: np.ndarray) -> float:
        """AngioEye formula: 1 - lambda2/lambda1."""
        s = np.asarray(s, dtype=float)
        if s.size < 2:
            return np.nan
        lam1 = float(s[0])
        lam2 = float(s[1])
        if not (np.isfinite(lam1) and np.isfinite(lam2) and lam1 > 0):
            return np.nan
        return float(1.0 - lam2 / lam1)

    # ------------------------------------------------------------------
    # Baseline / modal endpoints
    # ------------------------------------------------------------------

    def _compute_baseline_endpoints(
        self,
        mu: np.ndarray,
        x_full: np.ndarray,
        valid_column_mask: np.ndarray,
        T: np.ndarray,
        beat_period_valid: np.ndarray,
    ) -> dict:
        n_beats, n_branches, n_radii = mu.shape
        rms_x = np.full((n_beats, n_branches, n_radii), np.nan, dtype=float)
        rms_x[valid_column_mask] = np.sqrt(
            np.mean(x_full[:, valid_column_mask] ** 2, axis=0)
        )

        mean_pulsatile_ratio_bkr = np.full(
            (n_beats, n_branches, n_radii), np.nan, dtype=float
        )
        mean_pulsatile_ratio_bkr[valid_column_mask] = np.abs(
            mu[valid_column_mask]
        ) / (rms_x[valid_column_mask] + self.eps)

        tpr_b = self._median_kr_per_beat(rms_x, valid_column_mask)
        tpr = self._safe_nanmedian(tpr_b)
        mpr_b = self._median_kr_per_beat(mean_pulsatile_ratio_bkr, valid_column_mask)
        mpr = self._safe_nanmedian(mpr_b)
        abs_mu_b = self._median_kr_per_beat(np.abs(mu), valid_column_mask)
        abs_mu_acq = self._safe_nanmedian(abs_mu_b)
        mu_b = self._median_kr_per_beat(mu, valid_column_mask)
        mpr_prime = (
            float(abs_mu_acq / (tpr + self.eps))
            if np.isfinite(abs_mu_acq) and np.isfinite(tpr) and tpr > self.eps
            else np.nan
        )
        rho0_b = np.where(np.isfinite(tpr_b) & (tpr_b > self.eps), 1.0, np.nan)

        beatwise = {
            "mu_b": mu_b,
            "abs_mu_b": abs_mu_b,
            "TPR_b": tpr_b,
            "R0_b": tpr_b,  # EyeFlow alias
            "rho0_b": rho0_b,
            "mpr_b": mpr_b,
            "MPR_b": mpr_b,
            "spatial_mad_mu_b": self._spatial_mad_per_beat(mu, valid_column_mask),
        }
        acq = {
            "mu_acq": self._safe_nanmedian(mu_b),
            "beat_period_mean": self._safe_nanmean(T[0][beat_period_valid]),
            "beat_period_median": self._safe_nanmedian(T[0][beat_period_valid]),
            "beat_period_std": self._safe_nanstd(T[0][beat_period_valid]),
            "sigma_mu_beat": self._safe_nanstd(mu_b),
            "mad_mu_beat": self._safe_nanmad(mu_b),
            "TPR": tpr,
            "R0": tpr,  # EyeFlow alias of AngioEye TPR (two-stage median)
            "rho0": self._safe_nanmedian(rho0_b),
            "sigma_TPR_beat": self._safe_nanstd(tpr_b),
            "mad_TPR_beat": self._safe_nanmad(tpr_b),
            "sigma_R0_beat": self._safe_nanstd(tpr_b),
            "mad_R0_beat": self._safe_nanmad(tpr_b),
            "sigma_rho0_beat": self._safe_nanstd(rho0_b),
            "mad_rho0_beat": self._safe_nanmad(rho0_b),
            "cv_rho0_beat": self._safe_nancv(rho0_b),
            "mpr": mpr,
            "MPR": mpr,
            "sigma_mpr_beat": self._safe_nanstd(mpr_b),
            "mad_mpr_beat": self._safe_nanmad(mpr_b),
            "cv_mpr_beat": self._safe_nancv(mpr_b),
            "abs_mu_acq": abs_mu_acq,
            "mpr_prime": mpr_prime,
            "spatial_mad_mu_median_over_beats": self._safe_nanmedian(
                beatwise["spatial_mad_mu_b"]
            ),
            "spectrum_mode_count": 0,
        }
        return {
            "rms_x": rms_x,
            "mean_pulsatile_ratio_bkr": mean_pulsatile_ratio_bkr,
            "tpr_b": tpr_b,
            "tpr": tpr,
            "beatwise": beatwise,
            "acq": acq,
        }

    def _mode_component_rms(
        self, u: np.ndarray, scores: np.ndarray, valid_mask: np.ndarray
    ) -> np.ndarray:
        rms_u = float(np.sqrt(np.mean(np.asarray(u, dtype=float) ** 2)))
        comp = np.full(valid_mask.shape, np.nan, dtype=float)
        comp[valid_mask] = np.abs(np.asarray(scores, dtype=float)) * rms_u
        return comp

    @staticmethod
    def _mode_sign_should_flip(scores: np.ndarray) -> bool:
        med = LowRankWaveformDecompositionCalculator._safe_nanmedian(scores)
        return bool(np.isfinite(med) and med < 0)

    @staticmethod
    def _reconstruct_mode_sum(U_r: np.ndarray, scores_r: np.ndarray) -> np.ndarray:
        if U_r.size == 0 or scores_r.size == 0:
            return np.zeros((U_r.shape[0], scores_r.shape[1]), dtype=float)
        return U_r @ scores_r

    def _residual_t_bkr(
        self,
        x_full: np.ndarray,
        valid_column_mask: np.ndarray,
        residual_valid: np.ndarray,
    ) -> np.ndarray:
        residual = np.full_like(x_full, np.nan, dtype=float)
        residual[:, valid_column_mask] = residual_valid
        return residual

    def _compute_modal_endpoints(
        self,
        X: np.ndarray,
        x_full: np.ndarray,
        svd: dict,
        valid_column_mask: np.ndarray,
        tpr_b: np.ndarray,
        tpr: float,
        n_t: int,
        n_beats: int,
        n_branches: int,
        n_radii: int,
        beatwise: dict,
        acq: dict,
    ) -> dict:
        U = svd["U"]
        s = svd["s"]
        score_list = svd["score_list"]
        score_panel_bkr = svd["score_panel_bkr"]
        energy_fraction = svd["energy_fraction"]
        n_modes = svd["n_modes_panel"]

        rms_mode_panel = np.full_like(score_panel_bkr, np.nan, dtype=float)
        residual_rms_panel = np.full_like(score_panel_bkr, np.nan, dtype=float)
        residual_t_bkr_panel = np.full(
            (self.exported_modes, n_t, n_beats, n_branches, n_radii),
            np.nan,
            dtype=float,
        )

        for m in range(1, n_modes + 1):
            u_m = U[:, m - 1]
            scores_m = score_list[m - 1]
            rms_mode_panel[m - 1] = self._mode_component_rms(
                u=u_m, scores=scores_m, valid_mask=valid_column_mask
            )
            X_recon_m = self._reconstruct_mode_sum(
                U[:, :m], np.vstack(score_list[:m])
            )
            X_res_m = X - X_recon_m
            if m <= self.exported_modes:
                residual_t_bkr_panel[m - 1] = self._residual_t_bkr(
                    x_full=x_full,
                    valid_column_mask=valid_column_mask,
                    residual_valid=X_res_m,
                )
            residual_rms_bkr = np.full(
                (n_beats, n_branches, n_radii), np.nan, dtype=float
            )
            residual_rms_bkr[valid_column_mask] = np.sqrt(np.mean(X_res_m**2, axis=0))
            residual_rms_panel[m - 1] = residual_rms_bkr

            r_b = self._median_kr_per_beat(residual_rms_bkr, valid_column_mask)
            a_b = self._median_kr_per_beat(rms_mode_panel[m - 1], valid_column_mask)
            rho_b = np.where(
                np.isfinite(r_b) & np.isfinite(tpr_b) & (tpr_b > self.eps),
                r_b / (tpr_b + self.eps),
                np.nan,
            )
            R_m = self._safe_nanmedian(r_b)

            beatwise[f"A{m}_b"] = a_b
            beatwise[f"R{m}_b"] = r_b
            beatwise[f"rho{m}_b"] = rho_b
            beatwise[f"median_abs_a{m}_b"] = self._median_kr_per_beat(
                np.abs(score_panel_bkr[m - 1]), valid_column_mask
            )

            acq[f"A{m}"] = self._safe_nanmedian(a_b)
            acq[f"R{m}"] = R_m
            acq[f"rho{m}"] = (
                float(R_m / (tpr + self.eps))
                if np.isfinite(R_m) and np.isfinite(tpr) and tpr > self.eps
                else np.nan
            )
            acq[f"sigma_A{m}_beat"] = self._safe_nanstd(a_b)
            acq[f"mad_A{m}_beat"] = self._safe_nanmad(a_b)
            acq[f"cv_A{m}_beat"] = self._safe_nancv(a_b)
            acq[f"sigma_R{m}_beat"] = self._safe_nanstd(r_b)
            acq[f"mad_R{m}_beat"] = self._safe_nanmad(r_b)
            acq[f"cv_R{m}_beat"] = self._safe_nancv(r_b)
            acq[f"sigma_rho{m}_beat"] = self._safe_nanstd(rho_b)
            acq[f"mad_rho{m}_beat"] = self._safe_nanmad(rho_b)
            acq[f"cv_rho{m}_beat"] = self._safe_nancv(rho_b)
            acq[f"median_abs_a{m}"] = self._safe_nanmedian(
                beatwise[f"median_abs_a{m}_b"]
            )
            acq[f"spatial_mad_A{m}_median_over_beats"] = self._safe_nanmedian(
                self._spatial_mad_per_beat(rms_mode_panel[m - 1], valid_column_mask)
            )
            acq[f"spatial_mad_R{m}_median_over_beats"] = self._safe_nanmedian(
                self._spatial_mad_per_beat(residual_rms_bkr, valid_column_mask)
            )

        acq["eta1"] = float(energy_fraction[0]) if len(energy_fraction) >= 1 else np.nan
        acq["eta2"] = float(energy_fraction[1]) if len(energy_fraction) >= 2 else np.nan
        acq["eta12"] = (
            float(np.sum(energy_fraction[:2])) if len(energy_fraction) >= 1 else np.nan
        )
        acq["effective_rank"] = self._effective_rank(energy_fraction)
        acq["participation_ratio"] = self._participation_ratio(energy_fraction)
        acq["alpha"] = self._alpha(s, n_modes)
        acq["G1"] = self._g1(s)
        acq["spectrum_mode_count"] = int(np.asarray(s).size)

        return {
            "rms_mode_panel": rms_mode_panel,
            "residual_rms_panel": residual_rms_panel,
            "residual_t_bkr_panel": residual_t_bkr_panel,
        }

    def _run_joint_svd(
        self,
        X: np.ndarray,
        n_t: int,
        n_beats: int,
        n_branches: int,
        n_radii: int,
        n_valid_columns: int,
        valid_column_mask: np.ndarray,
    ) -> dict:
        U, s, Vt = np.linalg.svd(X, full_matrices=False)
        energy = s**2
        energy_fraction = energy / (np.sum(energy) + self.eps)
        n_modes = int(min(self.exported_modes, len(s)))

        score_list: list[np.ndarray] = []
        sign_flips = np.zeros((n_modes,), dtype=int)
        u_panel = np.full((n_t, self.exported_modes), np.nan, dtype=float)
        score_panel_flat = np.full(
            (self.exported_modes, n_valid_columns), np.nan, dtype=float
        )

        for m in range(n_modes):
            scores = s[m] * Vt[m, :]
            if self._mode_sign_should_flip(scores):
                U[:, m] *= -1.0
                Vt[m, :] *= -1.0
                scores *= -1.0
                sign_flips[m] = 1
            u_panel[:, m] = U[:, m]
            score_panel_flat[m, :] = scores
            score_list.append(scores)

        score_panel_bkr = np.full(
            (self.exported_modes, n_beats, n_branches, n_radii), np.nan, dtype=float
        )
        for m in range(n_modes):
            score_panel_bkr[m, valid_column_mask] = score_list[m]

        return {
            "U": U,
            "s": s,
            "Vt": Vt,
            "energy": energy,
            "energy_fraction": energy_fraction,
            "n_modes_panel": n_modes,
            "score_list": score_list,
            "sign_flips": sign_flips,
            "U_panel": u_panel,
            "score_panel_flat": score_panel_flat,
            "score_panel_bkr": score_panel_bkr,
        }

    def _svd_beat_panel(self, v_beat: np.ndarray, beat_period: float) -> dict:
        max_modes = self.exported_modes
        _, n_branches, n_radii = v_beat.shape
        mode_rms = np.full((max_modes, n_branches, n_radii), np.nan, dtype=float)
        residual_rms = np.full((max_modes, n_branches, n_radii), np.nan, dtype=float)
        total_rms = np.full((n_branches, n_radii), np.nan, dtype=float)
        mean_pulsatile_ratio = np.full((n_branches, n_radii), np.nan, dtype=float)
        effective_rank = np.nan
        participation_ratio = np.nan
        alpha = np.nan
        g1 = np.nan

        def _pack(valid_mask: np.ndarray) -> dict:
            return {
                "mode_rms": mode_rms,
                "residual_rms": residual_rms,
                "total_rms": total_rms,
                "mean_pulsatile_ratio": mean_pulsatile_ratio,
                "valid_mask": valid_mask,
                "effective_rank": effective_rank,
                "participation_ratio": participation_ratio,
                "alpha": alpha,
                "g1": g1,
                "singular_values": np.asarray([], dtype=float),
            }

        finite_fraction = np.mean(np.isfinite(v_beat), axis=0)
        valid_mask = finite_fraction >= float(self.min_valid_samples_fraction)
        beat_period_valid = np.isfinite(beat_period) and beat_period > 0
        if not beat_period_valid:
            valid_mask[:] = False
            return _pack(valid_mask)
        if int(np.sum(valid_mask)) < int(self.min_valid_columns):
            return _pack(valid_mask)

        mu, x_full = mean_subtract(v_beat)
        X = x_full[:, valid_mask]
        if X.size == 0:
            return _pack(valid_mask)

        total_rms[valid_mask] = np.sqrt(np.mean(X**2, axis=0))
        mean_pulsatile_ratio[valid_mask] = np.abs(mu[valid_mask]) / (
            total_rms[valid_mask] + self.eps
        )

        U, s, Vt = np.linalg.svd(X, full_matrices=False)
        energy = s**2
        energy_fraction = energy / (np.sum(energy) + self.eps)
        effective_rank = self._effective_rank(energy_fraction)
        participation_ratio = self._participation_ratio(energy_fraction)
        n_modes_alpha = int(min(self.exported_modes, len(s)))
        alpha = self._alpha(s, n_modes_alpha)
        g1 = self._g1(s)

        n_modes = int(min(max_modes, len(s)))
        scores_list: list[np.ndarray] = []
        for m in range(n_modes):
            scores = s[m] * Vt[m, :]
            if self._mode_sign_should_flip(scores):
                U[:, m] *= -1.0
                scores = -scores
            scores_list.append(scores)
            rms_u = float(np.sqrt(np.mean(U[:, m] ** 2)))
            mode_rms[m][valid_mask] = np.abs(scores) * rms_u
            X_recon = self._reconstruct_mode_sum(
                U[:, : m + 1], np.vstack(scores_list[: m + 1])
            )
            residual = X - X_recon
            residual_rms[m][valid_mask] = np.sqrt(np.mean(residual**2, axis=0))

        packed = _pack(valid_mask)
        packed["singular_values"] = np.asarray(s, dtype=float)
        return packed

    def per_beat_svd_panels(self, v_block: np.ndarray, T: np.ndarray) -> dict:
        T = normalize_periods(T)
        _, n_beats, n_branches, n_radii = v_block.shape
        max_modes = self.exported_modes
        mode_rms = np.full(
            (max_modes, n_beats, n_branches, n_radii), np.nan, dtype=float
        )
        residual_rms = np.full(
            (max_modes, n_beats, n_branches, n_radii), np.nan, dtype=float
        )
        total_rms = np.full((n_beats, n_branches, n_radii), np.nan, dtype=float)
        mean_pulsatile_ratio = np.full(
            (n_beats, n_branches, n_radii), np.nan, dtype=float
        )
        valid_mask = np.zeros((n_beats, n_branches, n_radii), dtype=bool)
        effective_rank_b = np.full((n_beats,), np.nan, dtype=float)
        participation_ratio_b = np.full((n_beats,), np.nan, dtype=float)
        alpha_b = np.full((n_beats,), np.nan, dtype=float)
        g1_b = np.full((n_beats,), np.nan, dtype=float)
        # Keep enough leading singular values for cohort spectrum plots.
        n_spectrum = max(int(max_modes), 12)
        singular_values_b = np.full((n_beats, n_spectrum), np.nan, dtype=float)

        for b in range(n_beats):
            beat = self._svd_beat_panel(v_block[:, b, :, :], beat_period=float(T[0, b]))
            mode_rms[:, b, :, :] = beat["mode_rms"]
            residual_rms[:, b, :, :] = beat["residual_rms"]
            total_rms[b, :, :] = beat["total_rms"]
            mean_pulsatile_ratio[b, :, :] = beat["mean_pulsatile_ratio"]
            valid_mask[b, :, :] = beat["valid_mask"]
            effective_rank_b[b] = beat["effective_rank"]
            participation_ratio_b[b] = beat["participation_ratio"]
            alpha_b[b] = beat["alpha"]
            g1_b[b] = beat["g1"]
            s_beat = np.asarray(beat.get("singular_values", []), dtype=float)
            n_keep = int(min(n_spectrum, s_beat.size))
            if n_keep > 0:
                singular_values_b[b, :n_keep] = s_beat[:n_keep]

        return {
            "mode_rms": mode_rms,
            "residual_rms": residual_rms,
            "total_rms": total_rms,
            "mean_pulsatile_ratio": mean_pulsatile_ratio,
            "valid_mask": valid_mask,
            "effective_rank_b": effective_rank_b,
            "participation_ratio_b": participation_ratio_b,
            "alpha_b": alpha_b,
            "g1_b": g1_b,
            "singular_values_b": singular_values_b,
        }

    def _compute_per_beat_endpoints(self, panels: dict) -> dict[str, np.ndarray]:
        valid_mask = panels["valid_mask"]
        tpr_b_pb = self._median_kr_per_beat(panels["total_rms"], valid_mask)
        mpr_b_pb = self._median_kr_per_beat(
            panels["mean_pulsatile_ratio"], valid_mask
        )
        a1_b_pb = self._median_kr_per_beat(panels["mode_rms"][0], valid_mask)
        a2_b_pb = self._median_kr_per_beat(panels["mode_rms"][1], valid_mask)
        r1_b_pb = self._median_kr_per_beat(panels["residual_rms"][0], valid_mask)
        r2_b_pb = self._median_kr_per_beat(panels["residual_rms"][1], valid_mask)
        rho1_b_pb = np.where(
            np.isfinite(r1_b_pb) & np.isfinite(tpr_b_pb) & (tpr_b_pb > self.eps),
            r1_b_pb / (tpr_b_pb + self.eps),
            np.nan,
        )
        rho2_b_pb = np.where(
            np.isfinite(r2_b_pb) & np.isfinite(tpr_b_pb) & (tpr_b_pb > self.eps),
            r2_b_pb / (tpr_b_pb + self.eps),
            np.nan,
        )
        rho0_b_pb = np.where(
            np.isfinite(tpr_b_pb) & (tpr_b_pb > self.eps), 1.0, np.nan
        )
        return {
            "A1_b_pb": a1_b_pb,
            "A2_b_pb": a2_b_pb,
            "R1_b_pb": r1_b_pb,
            "R2_b_pb": r2_b_pb,
            "rho1_b_pb": rho1_b_pb,
            "rho2_b_pb": rho2_b_pb,
            "TPR_b_pb": tpr_b_pb,
            "R0_b_pb": tpr_b_pb,
            "rho0_b_pb": rho0_b_pb,
            "MPR_b_pb": mpr_b_pb,
            "mpr_b_pb": mpr_b_pb,
            "effective_rank_b_pb": panels["effective_rank_b"],
            "participation_ratio_b_pb": panels["participation_ratio_b"],
            "alpha_b_pb": panels["alpha_b"],
            "G1_b_pb": panels["g1_b"],
            "singular_values_b": np.asarray(
                panels.get("singular_values_b", []), dtype=float
            ),
        }

    def _compute_representation(self, v_block: np.ndarray, T: np.ndarray) -> dict:
        T = normalize_periods(T)
        v_block = ensure_segment_shape(v_block, T)
        n_t, n_beats, n_branches, n_radii = v_block.shape
        if T.shape[1] != n_beats:
            raise ValueError(
                "Beat-period length mismatch: "
                f"T has {T.shape[1]} beats, waveform block has {n_beats} beats."
            )

        finite_fraction = np.mean(np.isfinite(v_block), axis=0)
        valid_column_mask = finite_fraction >= float(self.min_valid_samples_fraction)
        beat_period_valid = np.isfinite(T[0]) & (T[0] > 0)
        if np.any(~beat_period_valid):
            valid_column_mask &= beat_period_valid[:, None, None]

        n_total_columns = int(n_beats * n_branches * n_radii)
        n_valid_columns = int(np.sum(valid_column_mask))

        out: dict = {
            "shape": {
                "n_t": n_t,
                "n_beats": n_beats,
                "n_branches": n_branches,
                "n_radii": n_radii,
                "n_total_columns": n_total_columns,
                "n_valid_columns": n_valid_columns,
            },
            "valid_column_mask": valid_column_mask,
            "finite_fraction_per_column": finite_fraction,
            "beat_period_valid": beat_period_valid,
        }

        mu, x_full = mean_subtract(v_block)
        out["mu"] = mu
        out["x_full"] = x_full
        valid_counts_per_beat = np.sum(valid_column_mask, axis=(1, 2))
        out["valid_counts_per_beat"] = valid_counts_per_beat
        out["valid_fraction_per_beat"] = valid_counts_per_beat / float(
            max(1, n_branches * n_radii)
        )

        baseline = self._compute_baseline_endpoints(
            mu=mu,
            x_full=x_full,
            valid_column_mask=valid_column_mask,
            T=T,
            beat_period_valid=beat_period_valid,
        )
        out["rms_x"] = baseline["rms_x"]
        out["total_rms_bkr"] = baseline["rms_x"]
        out["mean_pulsatile_ratio_bkr"] = baseline["mean_pulsatile_ratio_bkr"]
        # EyeFlow packer compatibility alias
        out["mean_to_pulsatile_ratio_bkr"] = baseline["mean_pulsatile_ratio_bkr"]
        beatwise = baseline["beatwise"]
        acq = baseline["acq"]
        tpr_b = baseline["tpr_b"]
        tpr = baseline["tpr"]

        if n_valid_columns < int(self.min_valid_columns):
            out["beatwise"] = beatwise
            out["acq"] = acq
            out["svd_available"] = False
            out["svd_reason"] = "too_few_valid_columns"
            out["energy_fraction"] = np.asarray([], dtype=float)
            return out

        X = x_full[:, valid_column_mask]
        if X.size == 0:
            out["beatwise"] = beatwise
            out["acq"] = acq
            out["svd_available"] = False
            out["svd_reason"] = "empty_valid_matrix"
            out["energy_fraction"] = np.asarray([], dtype=float)
            return out

        svd = self._run_joint_svd(
            X=X,
            n_t=n_t,
            n_beats=n_beats,
            n_branches=n_branches,
            n_radii=n_radii,
            n_valid_columns=n_valid_columns,
            valid_column_mask=valid_column_mask,
        )
        out["svd_available"] = True
        out["svd_reason"] = "ok"
        out["X"] = X
        out["U"] = svd["U"]
        out["s"] = svd["s"]
        out["Vt"] = svd["Vt"]
        out["energy"] = svd["energy"]
        out["energy_fraction"] = svd["energy_fraction"]
        out["n_modes_panel"] = svd["n_modes_panel"]
        out["U_panel"] = svd["U_panel"]
        out["score_panel_flat"] = svd["score_panel_flat"]
        out["sign_flips"] = svd["sign_flips"]
        out["score_panel_bkr"] = svd["score_panel_bkr"]

        modal = self._compute_modal_endpoints(
            X=X,
            x_full=x_full,
            svd=svd,
            valid_column_mask=valid_column_mask,
            tpr_b=tpr_b,
            tpr=tpr,
            n_t=n_t,
            n_beats=n_beats,
            n_branches=n_branches,
            n_radii=n_radii,
            beatwise=beatwise,
            acq=acq,
        )
        out["rms_mode_panel"] = modal["rms_mode_panel"]
        out["residual_rms_panel"] = modal["residual_rms_panel"]
        out["residual_t_bkr_panel"] = modal["residual_t_bkr_panel"]
        out["beatwise"] = beatwise
        out["acq"] = acq
        return out


__all__ = [
    "LowRankWaveformDecompositionCalculator",
    "ensure_segment_shape",
    "mean_subtract",
    "normalize_periods",
]
