from epicirc.geneset.lock import build_lock


def test_build_lock_hashes_real_config():
    lock = build_lock("config/gene_sets")
    assert "circadian_core" in lock["gene_sets"]
    core = lock["gene_sets"]["circadian_core"]
    assert core["n_genes"] >= 20
    assert len(core["sha256"]) == 64
    # canonical clock genes present
    assert "CLOCK" in core["symbols"] and "ARNTL" in core["symbols"]


def test_lock_audits_overlap_with_positive_control():
    lock = build_lock("config/gene_sets")
    # circadian core and the epilepsy positive control are disjoint by design -> no overlap key,
    # or if any gene were shared it must be surfaced transparently.
    for key, shared in lock["overlaps"].items():
        assert key.startswith("circadian_core∩")
        assert isinstance(shared, list) and shared
