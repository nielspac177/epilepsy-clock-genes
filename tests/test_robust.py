import pytest

from epicirc.mr.robust import mr_presso_lite, steiger_filter, weighted_median


def _clean(true_beta=0.5):
    bx = [0.10, 0.20, 0.15, 0.30, 0.25, 0.18, 0.22, 0.12, 0.28, 0.16]
    byse = [0.02] * len(bx)
    by = [true_beta * x for x in bx]
    return bx, by, byse


def test_weighted_median_recovers_true_effect():
    bx, by, byse = _clean(0.5)
    est = weighted_median(bx, by, byse, n_boot=200)
    assert est.beta == pytest.approx(0.5, abs=0.02)


def test_presso_recovers_and_flags_outlier():
    bx, by, byse = _clean(0.5)
    by[3] += 0.5  # inject a strong outlier at instrument 3
    est = mr_presso_lite(bx, by, byse)
    assert est.extra["n_outliers"] >= 1
    assert est.beta == pytest.approx(0.5, abs=0.05)  # corrected estimate ignores the outlier


def test_steiger_keeps_valid_direction():
    bx, by, byse = _clean(0.4)
    eaf = [0.3] * len(bx)
    est = steiger_filter(bx, byse, eaf, 400000, by, byse, 50000)
    # exposure effects are larger than outcome effects (beta<1) -> all kept, correct direction
    assert est.extra["n_removed"] == 0
    assert est.beta == pytest.approx(0.4, abs=0.02)
