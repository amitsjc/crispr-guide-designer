"""Tests for the gRNA designer."""

from crispr_designer.guide_designer import (
    design_guides, _gc_fraction, _homopolymer_run, _quality_flags,
    GC_MIN, GC_MAX, GuideRNA
)
from crispr_designer.pam_finder import GUIDE_LENGTH


def _make_seq(spacer: str, pam: str = "AGG", padding: str = "A") -> str:
    return spacer + pam + padding * 10


# --- helpers ---

def test_gc_fraction_pure_gc():
    assert _gc_fraction("GGCC") == 1.0

def test_gc_fraction_pure_at():
    assert _gc_fraction("AATT") == 0.0

def test_gc_fraction_mixed():
    assert _gc_fraction("ATGC") == 0.5

def test_gc_fraction_empty():
    assert _gc_fraction("") == 0.0

def test_homopolymer_run_single_char():
    assert _homopolymer_run("AAAA") == 4

def test_homopolymer_run_mixed():
    assert _homopolymer_run("AAGGG") == 3

def test_homopolymer_run_no_run():
    assert _homopolymer_run("ATGC") == 1

def test_quality_flags_poly_t():
    assert any("poly-T" in f for f in _quality_flags("ATTTTATTTTATCGATCGAT"))

def test_quality_flags_homopolymer():
    flags = _quality_flags("AAAAAAATCGATCGATCGAT")
    assert any("homopolymer" in f for f in flags)

def test_quality_flags_clean():
    assert _quality_flags("ATGCATGCATGCATGCATGC") == []


# --- design_guides ---

def test_gc_filter_excludes_low_gc():
    # All-AT spacer: GC=0%, should be excluded
    spacer = "A" * GUIDE_LENGTH
    seq = _make_seq(spacer)
    guides = design_guides(seq, gc_min=0.4, gc_max=0.7)
    assert all(g.gc_content >= 0.4 for g in guides)

def test_gc_filter_excludes_high_gc():
    # All-GC spacer: GC=100%, should be excluded
    spacer = "G" * GUIDE_LENGTH
    seq = _make_seq(spacer)
    guides = design_guides(seq, gc_min=0.4, gc_max=0.7)
    assert all(g.gc_content <= 0.7 for g in guides)

def test_good_guide_passes():
    # 50% GC spacer should pass default filters
    spacer = "ATGCATGCATGCATGCATGC"  # exactly 50% GC, 20nt
    seq = _make_seq(spacer)
    guides = design_guides(seq, gc_min=0.4, gc_max=0.7)
    assert any(g.spacer == spacer for g in guides)

def test_exclude_flagged_removes_homopolymer():
    spacer = "AAAAAATGCATGCATGCATG"  # 6 A's in a row
    seq = _make_seq(spacer)
    guides_all = design_guides(seq, exclude_flagged=False)
    guides_clean = design_guides(seq, exclude_flagged=True)
    flagged = [g for g in guides_all if g.flags]
    # All flagged guides should be absent from clean list
    for g in flagged:
        assert g.spacer not in [c.spacer for c in guides_clean]

def test_guide_positions_sorted():
    spacer = "ATGCATGCATGCATGCATGC"
    seq = _make_seq(spacer) * 3
    guides = design_guides(seq)
    positions = [g.position for g in guides]
    assert positions == sorted(positions)

def test_vegfa_produces_guides():
    with open("data/vegfa.fasta") as f:
        seq = "".join(l.strip() for l in f if not l.startswith(">"))
    guides = design_guides(seq)
    assert len(guides) > 0
    assert all(0.4 <= g.gc_content <= 0.7 for g in guides)
