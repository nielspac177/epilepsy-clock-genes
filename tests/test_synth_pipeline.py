"""Regression tests for the synthetic end-to-end path (no external tools)."""
from epicirc import synth
from epicirc.io import read_sumstats, write_sumstats
from epicirc.magma.run_gene import aggregate_genes
from epicirc.magma.run_set import competitive_test


def test_synth_is_deterministic():
    a = synth.generate("focal_epilepsy")
    b = synth.generate("focal_epilepsy")
    assert [r.beta for r in a] == [r.beta for r in b]
    assert a[0].snp == "rs000000"


def test_ground_truth_gge_more_circadian_than_focal(tmp_path):
    # GGE has a larger injected circadian delta than focal -> larger circadian gene-set beta.
    def core_beta(name):
        p = tmp_path / f"{name}.sumstats.gz"
        write_sumstats(synth.generate(name), p)
        genes = [(g, z, c) for g, _chr, _n, z, c in aggregate_genes(str(p))]
        in_set = {g for g, _z, c in genes if c == 1}
        return competitive_test(genes, in_set)[1]  # beta

    assert core_beta("generalized_gge") > core_beta("focal_epilepsy")


def test_circadian_column_roundtrips(tmp_path):
    p = tmp_path / "x.sumstats.gz"
    write_sumstats(synth.generate("chronotype_morning"), p)
    rows = list(read_sumstats(p))
    assert any(r["CIRCADIAN"] == "1" for r in rows)
    assert {"SNP", "BETA", "SE", "GENE"} <= set(rows[0])
