from epicirc.mr.coloc import coloc


def _null_snp():
    return {"e_p": 0.6, "g_beta": 0.002, "g_se": 0.01}


def _build(strong_eqtl_at, strong_gwas_at, n=25):
    eqtl, gwas = {}, {}
    for i in range(n):
        k = ("1", 1000 + i)
        e_p = 1e-30 if i == strong_eqtl_at else 0.6
        g_beta = 0.20 if i == strong_gwas_at else 0.002
        eqtl[k] = e_p
        gwas[k] = (g_beta, 0.01, 0.3)
    return eqtl, gwas


def test_shared_variant_gives_high_h4():
    eqtl, gwas = _build(strong_eqtl_at=10, strong_gwas_at=10)
    r = coloc(eqtl, gwas)
    assert r["PP.H4"] > 0.9


def test_distinct_variants_give_high_h3():
    eqtl, gwas = _build(strong_eqtl_at=5, strong_gwas_at=18)
    r = coloc(eqtl, gwas)
    assert r["PP.H3"] > 0.5
    assert r["PP.H4"] < 0.2


def test_too_few_snps_returns_nan():
    import math
    eqtl = {("1", i): 0.5 for i in range(3)}
    gwas = {("1", i): (0.01, 0.01, 0.3) for i in range(3)}
    assert math.isnan(coloc(eqtl, gwas)["PP.H4"])
