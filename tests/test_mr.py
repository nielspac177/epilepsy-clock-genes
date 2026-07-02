import pytest

from epicirc.mr.two_sample import ivw, mr_egger


def _instruments(true_beta, pleiotropy=0.0):
    bx = [0.10, 0.20, 0.15, 0.30, 0.25, 0.18, 0.22]
    sey = [0.02, 0.02, 0.03, 0.02, 0.025, 0.02, 0.03]
    by = [true_beta * x + pleiotropy for x in bx]
    return bx, by, sey


def test_ivw_recovers_true_effect():
    bx, by, sey = _instruments(0.7)
    est = ivw(bx, by, sey)
    assert est.beta == pytest.approx(0.7, rel=1e-6)
    assert est.method == "ivw"
    assert est.extra["cochran_q"] == pytest.approx(0.0, abs=1e-9)  # no heterogeneity


def test_egger_slope_matches_and_intercept_zero_without_pleiotropy():
    bx, by, sey = _instruments(0.7)
    est = mr_egger(bx, by, sey)
    assert est.beta == pytest.approx(0.7, rel=1e-6)
    assert est.extra["intercept"] == pytest.approx(0.0, abs=1e-9)


def test_egger_intercept_detects_directional_pleiotropy():
    bx, by, sey = _instruments(0.7, pleiotropy=0.05)
    est = mr_egger(bx, by, sey)
    assert est.extra["intercept"] == pytest.approx(0.05, abs=1e-6)
    # slope still recovers the causal effect once pleiotropy is absorbed by the intercept
    assert est.beta == pytest.approx(0.7, rel=1e-6)


def test_mr_input_validation():
    with pytest.raises(ValueError):
        ivw([0.1], [0.1], [0.1])                    # too few instruments
    with pytest.raises(ValueError):
        ivw([0.1, 0.2], [0.1, 0.2], [0.0, 0.1])     # non-positive SE
    with pytest.raises(ValueError):
        mr_egger([0.1, 0.2], [0.1, 0.2], [0.1, 0.1])  # egger needs >=3
