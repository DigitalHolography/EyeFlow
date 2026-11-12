"""Population STD across existing beats after phase-window median reduction."""

import warnings

import numpy as np

from pipelines.low_freq_spectra.per_beat import beat_spectrum_std
from calculations.math.phase_windows import early_late_phase_windows


def _raw():
    values = np.full((4, 3, 3), np.nan, dtype=np.float32)
    values[:, 0, :] = [[0, 2, 4], [2, 6, 10], [10, 18, 26], [14, 22, 30]]
    values[:, 1, 0] = [1, np.nan, 3, np.nan]
    return values


def test_std_uses_per_beat_phase_medians_population_weighting_and_nan_counts():
    raw = _raw()
    original = raw.copy()
    windows = {"whole": (np.ones(4, bool), [0, 1]), **early_late_phase_windows(np.arange(4) / 4)}
    for name, expected in (("whole", np.sqrt(24)), ("early", np.sqrt(6)), ("late", np.sqrt(128 / 3))):
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            std, counts = beat_spectrum_std(raw, windows[name][0])
        np.testing.assert_allclose(std[:2], [expected, 0])
        assert np.isnan(std[2])
        np.testing.assert_array_equal(counts, [3, 1, 0])
        assert std.shape == (3,) and std.dtype == np.float32
    np.testing.assert_equal(raw, original)
