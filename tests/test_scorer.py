"""Tests for the on-target activity scorer."""

from crispr_designer.scorer import score_guide, rank_guides, _gc_score, _position_score
from crispr_designer.guide_designer import GuideRNA, design_guides
from crispr_designer.pam_finder import PAMSite


def _make_guide(spacer: str, gc: float = 0.50, flags=None) -> GuideRNA:
    site = PAMSite(position=20, strand="+", pam_seq="AGG", protospacer=spacer)
    return GuideRNA(spacer=spacer, pam_site=site, gc_content=gc, flags=flags or [])


# --- _gc_score ---

def test_gc_score_optimal():
    # 55% GC = optimal, should give maximum bonus
    assert _gc_score(0.55) > _gc_score(0.40)

def test_gc_score_extreme_low():
    # 40% GC (boundary) → some score
    s40 = _gc_score(0.40)
    s55 = _gc_score(0.55)
    assert s40 < s55

def test_gc_score_never_negative():
    assert _gc_score(0.0) >= 0.0
    assert _gc_score(1.0) >= 0.0


# --- _position_score ---

def test_position_score_g_at_20_boosts():
    # G at last position (index 19) should boost score
    spacer_with_g  = "ATCGATCGATCGATCGATCG"   # ends with G
    spacer_with_t  = "ATCGATCGATCGATCGATCT"   # ends with T
    assert _position_score(spacer_with_g) > _position_score(spacer_with_t)

def test_position_score_a_run_penalizes():
    # A at positions 2-4 should reduce score
    spacer_bad  = "GCAAATCGATCGATCGATCG"  # AAA at positions 2-4
    spacer_good = "GCGCATCGATCGATCGATCG"
    assert _position_score(spacer_bad) < _position_score(spacer_good)


# --- score_guide ---

def test_score_in_range():
    guide = _make_guide("ATGCATGCATGCATGCATGC", gc=0.50)
    s = score_guide(guide)
    assert 0.0 <= s <= 1.0

def test_flagged_guide_scores_lower():
    clean = _make_guide("ATGCATGCATGCATGCATGC", gc=0.50, flags=[])
    flagged = _make_guide("ATGCATGCATGCATGCATGC", gc=0.50, flags=["homopolymer run of 5"])
    assert score_guide(flagged) < score_guide(clean)

def test_optimal_gc_scores_higher():
    guide_50 = _make_guide("ATGCATGCATGCATGCATGC", gc=0.50)
    guide_40 = _make_guide("ATGCATGCATGCATGCATGC", gc=0.40)
    assert score_guide(guide_50) >= score_guide(guide_40)


# --- rank_guides ---

def test_rank_guides_sorted_descending():
    with open("data/vegfa.fasta") as f:
        seq = "".join(l.strip() for l in f if not l.startswith(">"))
    guides = design_guides(seq)
    ranked = rank_guides(guides)
    scores = [s for _, s in ranked]
    assert scores == sorted(scores, reverse=True)

def test_rank_guides_returns_all():
    with open("data/vegfa.fasta") as f:
        seq = "".join(l.strip() for l in f if not l.startswith(">"))
    guides = design_guides(seq)
    ranked = rank_guides(guides)
    assert len(ranked) == len(guides)

def test_rank_guides_empty():
    assert rank_guides([]) == []
