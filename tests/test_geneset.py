import textwrap

import pytest

from epicirc.geneset.curate import (
    load_geneset,
    matched_negative_controls,
    overlap,
    validate_against_annotation,
)

CORE_TSV = textwrap.dedent("""\
    # comment line ignored
    symbol\tentrez_id\trationale
    ARNTL\t406\tpositive arm
    CLOCK\t9575\tpositive arm
    PER1\t5187\tnegative arm
""")


@pytest.fixture
def core_file(tmp_path):
    p = tmp_path / "core.tsv"
    p.write_text(CORE_TSV)
    return p


def test_load_and_hash_is_order_invariant(tmp_path, core_file):
    gs = load_geneset(core_file)
    assert set(gs.symbols) == {"ARNTL", "CLOCK", "PER1"}
    # reordering rows must not change the hash (canonicalized by sorting)
    reordered = tmp_path / "reordered.tsv"
    reordered.write_text("symbol\tx\nPER1\t1\nARNTL\t2\nCLOCK\t3\n")
    assert load_geneset(reordered).sha256 == gs.sha256


def test_duplicate_symbols_rejected(tmp_path):
    p = tmp_path / "dup.tsv"
    p.write_text("symbol\nARNTL\nARNTL\n")
    with pytest.raises(ValueError):
        load_geneset(p)


def test_validate_against_annotation_flags_unknown(tmp_path, core_file):
    gs = load_geneset(core_file)
    gene_loc = tmp_path / "genes.loc"
    # MAGMA layout: ENTREZ CHR START STOP STRAND SYMBOL  (PER1 intentionally missing)
    gene_loc.write_text("406 11 13299273 13408812 + ARNTL\n9575 4 55435158 55547958 + CLOCK\n")
    missing = validate_against_annotation(gs, gene_loc)
    assert missing == ["PER1"]


def test_overlap_reports_shared_genes(tmp_path, core_file):
    gs = load_geneset(core_file)
    other = tmp_path / "epi.tsv"
    other.write_text("symbol\nCLOCK\nSCN1A\n")
    assert overlap(gs, load_geneset(other)) == ("CLOCK",)


def test_matched_controls_are_deterministic_and_sized():
    gs = load_geneset  # noqa: F841 (import check)
    from epicirc.geneset.curate import GeneSet

    target = GeneSet("t", ("ARNTL", "CLOCK", "PER1"))
    universe = {f"G{i}": 10000 + i * 100 for i in range(200)}
    universe.update({"ARNTL": 10500, "CLOCK": 10600, "PER1": 10700})
    a = matched_negative_controls(target, universe, n_sets=5, seed=42)
    b = matched_negative_controls(target, universe, n_sets=5, seed=42)
    assert a == b                              # deterministic given seed
    assert len(a) == 5 and all(len(s) == 3 for s in a)
    assert all(g not in target.symbols for s in a for g in s)  # excludes target genes
