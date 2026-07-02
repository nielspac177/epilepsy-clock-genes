from epicirc.mr.clump import Variant, distance_clump


def test_keeps_lead_and_drops_neighbours():
    v = [
        Variant("rs1", "1", 1_000_000, 1e-20),   # lead
        Variant("rs2", "1", 1_500_000, 1e-9),    # within 1Mb of rs1 -> dropped
        Variant("rs3", "1", 3_000_000, 1e-12),   # >1Mb away -> kept
        Variant("rs4", "2", 1_000_000, 1e-8),    # other chr -> kept
    ]
    kept = {x.snp for x in distance_clump(v, window_kb=1000)}
    assert kept == {"rs1", "rs3", "rs4"}


def test_strongest_wins_within_window():
    v = [
        Variant("weak", "5", 500_000, 1e-8),
        Variant("strong", "5", 600_000, 1e-30),  # 100kb away, stronger
    ]
    kept = distance_clump(v, window_kb=1000)
    assert len(kept) == 1 and kept[0].snp == "strong"


def test_empty():
    assert distance_clump([], window_kb=1000) == []
