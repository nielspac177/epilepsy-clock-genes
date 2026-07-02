import math

import pytest

from epicirc.stats.heterogeneity import heterogeneity_test
from epicirc.stats.multipletesting import benjamini_hochberg, bonferroni


def test_identical_effects_give_null_heterogeneity():
    r = heterogeneity_test(0.3, 0.1, 0.3, 0.1, rho=0.0)
    assert r.delta == 0.0
    assert r.p_two_sided == pytest.approx(1.0)


def test_shared_control_correlation_recovers_power():
    # For a case-case DIFFERENCE with shared controls, positive rho shrinks the difference
    # variance (control freq cancels), so accounting for rho recovers power vs the rho=0 default.
    indep = heterogeneity_test(0.5, 0.1, 0.2, 0.1, rho=0.0)
    shared = heterogeneity_test(0.5, 0.1, 0.2, 0.1, rho=0.6)
    assert shared.se_delta < indep.se_delta
    assert shared.p_two_sided < indep.p_two_sided


def test_missing_rho_warns_and_defaults_to_zero():
    with pytest.warns(UserWarning):
        r = heterogeneity_test(0.5, 0.1, 0.2, 0.1)
    assert r.rho_assumed is True
    assert r.rho_used == 0.0


def test_heterogeneity_rejects_bad_inputs():
    with pytest.raises(ValueError):
        heterogeneity_test(0.1, 0.0, 0.1, 0.1, rho=0.0)     # zero SE
    with pytest.raises(ValueError):
        heterogeneity_test(0.1, 0.1, 0.1, 0.1, rho=1.5)     # rho out of range


def test_known_z_value():
    # delta=0.3, se=sqrt(0.02) -> z ~ 2.121
    r = heterogeneity_test(0.5, 0.1, 0.2, 0.1, rho=0.0)
    assert r.z == pytest.approx(0.3 / math.sqrt(0.02), rel=1e-9)


def test_bonferroni_scales_and_caps():
    assert bonferroni([0.01, 0.5]) == [0.02, 1.0]
    assert bonferroni([]) == []


def test_bh_monotone_and_ordered():
    p = [0.001, 0.008, 0.039, 0.041, 0.9]
    q = benjamini_hochberg(p)
    assert all(0.0 <= x <= 1.0 for x in q)
    # q-values non-decreasing along the sorted p ranking
    order = sorted(range(len(p)), key=lambda i: p[i])
    sorted_q = [q[i] for i in order]
    assert sorted_q == sorted(sorted_q)
    # smallest p gets q = p * m / 1
    assert q[0] == pytest.approx(min(0.001 * 5 / 1, 1.0))


def test_bh_rejects_out_of_range():
    with pytest.raises(ValueError):
        benjamini_hochberg([0.1, 1.2])
