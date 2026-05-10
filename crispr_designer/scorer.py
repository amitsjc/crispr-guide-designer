"""On-target activity scorer — ranks guide RNAs by predicted cutting efficiency."""

from crispr_designer.guide_designer import GuideRNA

# Position-specific nucleotide weights derived from Doench et al. 2014
# (Rule Set 1 — the first widely adopted CRISPR scoring model).
# Keys: (position_0based, nucleotide) → weight adjustment.
# Positive = promotes activity, negative = reduces activity.
POSITION_WEIGHTS: dict[tuple[int, str], float] = {
    # Strong preference for G at positions 17-20 (seed region, 3' end)
    (16, "G"): 0.10, (17, "G"): 0.15, (18, "G"): 0.20, (19, "G"): 0.25,
    # Prefer C at position 1 (5' end of spacer)
    (0,  "C"): 0.05,
    # Avoid A at positions 3-5
    (2,  "A"): -0.10, (3, "A"): -0.08, (4, "A"): -0.06,
    # Avoid T at positions 16-18 (seed region)
    (15, "T"): -0.10, (16, "T"): -0.15, (17, "T"): -0.15,
    # Prefer A or T at the cut site region (positions 10-12)
    (9,  "A"): 0.05, (10, "A"): 0.05, (11, "T"): 0.05,
}


def _position_score(spacer: str) -> float:
    """Sum position-specific weights for a 20nt spacer."""
    spacer = spacer.upper()
    score = 0.0
    for (pos, nuc), weight in POSITION_WEIGHTS.items():
        if pos < len(spacer) and spacer[pos] == nuc:
            score += weight
    return score


def _gc_score(gc_content: float) -> float:
    """
    Reward GC content near the optimal 50-60% range.
    Penalize guides far from this window. Based on the observation
    that mid-range GC guides show highest on-target efficiency.
    """
    optimal = 0.55
    deviation = abs(gc_content - optimal)
    return max(0.0, 0.20 - deviation)   # max bonus 0.20, decays linearly


def score_guide(guide: GuideRNA) -> float:
    """
    Compute an on-target activity score for a guide RNA (0.0 – 1.0).

    Combines:
      - Position-specific nucleotide weights (Doench 2014 rules)
      - GC content reward near optimal 50-60%
      - Flag penalty for quality issues

    Higher score = predicted higher on-target cutting efficiency.
    """
    base = 0.50   # baseline for any GC-passing guide
    pos  = _position_score(guide.spacer)
    gc   = _gc_score(guide.gc_content)
    flag_penalty = len(guide.flags) * 0.10

    raw = base + pos + gc - flag_penalty
    return round(max(0.0, min(1.0, raw)), 4)


def rank_guides(guides: list[GuideRNA]) -> list[tuple[GuideRNA, float]]:
    """
    Score all guides and return them sorted best-first.
    Returns list of (guide, score) tuples.
    """
    scored = [(g, score_guide(g)) for g in guides]
    scored.sort(key=lambda x: x[1], reverse=True)
    return scored
