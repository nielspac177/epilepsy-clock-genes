"""Correctness tests for the two-GWAS pairwise coloc and the 3-trait moloc enumeration."""
import math

from epicirc.mr.coloc_multi import coloc_pairwise, moloc3, _canonical_configs, _config_name


def _region(strong_at, beta_strong=0.30, se=0.02, n=30):
    """A GWAS region: strong beta at index `strong_at` (or None), null elsewhere."""
    d = {}
    for i in range(n):
        k = ("17", 8_040_000 + i * 500)
        beta = beta_strong if i == strong_at else 0.001
        d[k] = (beta, se)
    return d


def test_15_configurations():
    cfgs = _canonical_configs()
    assert len(cfgs) == 15
    names = {_config_name(c) for c in cfgs}
    # sanity: the fully-shared and fully-independent configs are present
    assert "abc" in names
    assert "H0" in names
    assert "a_b_c" in names  # all three distinct singletons


def test_pairwise_shared_variant_high_h4():
    a = _region(10)          # strong at same SNP
    b = _region(10)
    r = coloc_pairwise(a, b)
    assert r["PP.H4"] > 0.9, r


def test_pairwise_distinct_variants_high_h3():
    a = _region(5)
    b = _region(20)          # strong at a different SNP
    r = coloc_pairwise(a, b)
    assert r["PP.H3"] > 0.7, r
    assert r["PP.H4"] < 0.1, r


def test_pairwise_one_null_high_h1_or_h2():
    a = _region(10)
    b = _region(None)        # b has no signal
    r = coloc_pairwise(a, b)
    assert r["PP.H1"] > 0.5, r   # only trait a associated


def test_moloc3_all_shared():
    a, b, c = _region(12), _region(12), _region(12)
    r = moloc3(a, b, c)
    assert r["PP_abc"] > 0.8, r
    assert r["PP_none_shared"] < 0.2, r


def test_moloc3_all_distinct():
    a, b, c = _region(3), _region(13), _region(25)
    r = moloc3(a, b, c)
    assert r["PP_abc"] < 0.05, r
    # the fully-distinct configuration should dominate the shared ones
    assert r["cfg.a_b_c"] > r["cfg.abc"], r


def test_moloc3_two_share_one_distinct():
    a, b, c = _region(12), _region(12), _region(25)   # a,b share; c separate
    r = moloc3(a, b, c)
    # posterior mass should favor the (ab)c configuration over full abc and full-distinct
    assert r["cfg.ab_c"] > r["cfg.abc"], r
    assert r["cfg.ab_c"] > r["cfg.a_b_c"], r


def test_posteriors_sum_to_one():
    a, b, c = _region(12), _region(12), _region(12)
    r = moloc3(a, b, c)
    total = sum(v for k, v in r.items() if k.startswith("cfg."))
    assert abs(total - 1.0) < 1e-6, total
