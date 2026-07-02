import pytest

from epicirc.data.harmonize import harmonize, is_palindromic


def test_direct_match_keeps_sign():
    h = harmonize(0.4, "A", "G", ref_effect="A", ref_other="G")
    assert h.action == "match" and h.beta == 0.4 and h.usable


def test_allele_swap_flips_sign():
    h = harmonize(0.4, "G", "A", ref_effect="A", ref_other="G")
    assert h.action == "flip_sign" and h.beta == -0.4 and h.usable


def test_strand_flip_recognized():
    # effect C/T against reference G/A (complements)
    h = harmonize(0.25, "C", "T", ref_effect="G", ref_other="A")
    assert h.action == "strand_flip" and h.beta == 0.25 and h.usable


def test_strand_flip_with_swap():
    h = harmonize(0.25, "T", "C", ref_effect="G", ref_other="A")
    assert h.action == "strand_flip_and_sign" and h.beta == -0.25 and h.usable


def test_palindromic_dropped_by_default():
    assert is_palindromic("A", "T")
    h = harmonize(0.3, "A", "T", ref_effect="A", ref_other="T")
    assert h.action == "ambiguous" and not h.usable


def test_incompatible_alleles_flagged_unusable():
    h = harmonize(0.3, "A", "C", ref_effect="A", ref_other="G")
    assert h.action == "incompatible" and not h.usable
