import textwrap

from epicirc.analysis.geneset_enrich import gene_top_chi2, load_genes, matched_null


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(textwrap.dedent(text))
    return str(p)


def test_gene_top_chi2_picks_max_in_window(tmp_path):
    genes = _write(tmp_path, "genes.tsv", """\
        GENEA\t1\t1000\t2000
        GENEB\t1\t100000\t101000
    """)
    # ILAE .tbl: cols CHR BP SNP A1 A2 FRQ FreqSE Neff Z P Dir Beta SE
    tbl = _write(tmp_path, "x.tbl", """\
        CHR\tBP\tSNP\tA1\tA2\tFRQ\tFSE\tN\tZ\tP\tD\tB\tSE
        1\t1500\trs1\ta\tg\t0.3\t0\t1000\t2.0\t0.05\t+\t0.1\t0.05
        1\t1600\trs2\ta\tg\t0.3\t0\t1000\t4.0\t1e-4\t+\t0.2\t0.05
        1\t500000\trs3\ta\tg\t0.3\t0\t1000\t9.0\t1e-9\t+\t0.3\t0.05
    """)
    by_chr, length = load_genes(genes, flank_kb=1)  # ±1kb
    top, nsnp = gene_top_chi2(tbl, by_chr)
    assert top["GENEA"] == 16.0        # max(2^2, 4^2) within window
    assert nsnp["GENEA"] == 2
    assert "GENEB" not in top          # rs3 at 500k is outside GENEB's window


def test_matched_null_is_deterministic_and_returns_pvalue():
    top = {f"G{i}": float(i % 7) for i in range(400)}
    length = {f"G{i}": 1000 + i for i in range(400)}
    nsnp = {f"G{i}": 5 + (i % 20) for i in range(400)}
    clock = [f"G{i}" for i in range(0, 60, 3)]
    a = matched_null(top, length, nsnp, clock, n_null=500, seed=1)
    b = matched_null(top, length, nsnp, clock, n_null=500, seed=1)
    assert a == b
    assert 0 < a["emp_p"] <= 1
    assert a["obs_lo"] <= a["obs"] <= a["obs_hi"]      # bootstrap CI brackets the estimate
    assert a["null_lo"] <= a["null_mean"] <= a["null_hi"]
