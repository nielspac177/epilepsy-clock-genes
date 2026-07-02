import pytest

from epicirc.stats.prereg import effective_n, expected_h2_zscore, plan_cells


def test_effective_n_symmetric_and_bounded():
    # balanced design -> N_eff = N per group; huge controls -> approaches 4*Ncase
    assert effective_n(1000, 1000) == pytest.approx(2000.0)
    assert effective_n(1000, 10**9) == pytest.approx(4000.0, rel=1e-3)


def test_effective_n_rejects_nonpositive():
    with pytest.raises(ValueError):
        effective_n(0, 100)


def test_expected_z_orders_by_power_like_the_review():
    # GGE and focal should clear the gate; JAE/CAE should not (matches adversarial power analysis).
    gge = expected_h2_zscore(effective_n(7407, 52538), 0.32)
    focal = expected_h2_zscore(effective_n(16384, 52538), 0.18)
    cae = expected_h2_zscore(effective_n(1310, 52538), 0.40)
    jae = expected_h2_zscore(effective_n(607, 52538), 0.40)
    assert focal > 4 and gge > 4
    assert cae < 4 and jae < 4
    assert jae < cae < gge  # power ordering


def test_plan_cells_flags_estimable(tmp_path):
    manifest = {
        "epilepsy": {
            "generalized_gge": {"n_cases": 7407, "n_controls": 52538},
            "jae": {"n_cases": 607, "n_controls": 52538},
        }
    }
    cells = {c.phenotype: c for c in plan_cells(manifest)}
    assert cells["generalized_gge"].estimable is True
    assert cells["jae"].estimable is False
