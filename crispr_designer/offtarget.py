"""Off-target risk estimator — scores guide RNAs for potential unintended cuts."""

from dataclasses import dataclass
from crispr_designer.guide_designer import GuideRNA
from crispr_designer.pam_finder import find_pam_sites

# The seed region is the 12nt immediately adjacent to the PAM (positions 8-20).
# Mismatches here are much better tolerated than in the seed — a guide with
# even 1-2 mismatches in the seed can still cut off-target sites.
SEED_START = 8      # 0-based, counting from 5' end of the 20nt spacer
SEED_END   = 20

# Mismatch penalty weights by position within seed (heavier near PAM)
# Position 0 = most 5' of spacer, position 19 = adjacent to PAM
MISMATCH_WEIGHTS = {
    i: 1.0 + (i - SEED_START) * 0.15
    for i in range(SEED_START, SEED_END)
}


@dataclass
class OffTargetHit:
    """A potential off-target cut site in a background sequence."""
    sequence_id: str
    position: int
    strand: str
    candidate: str      # the off-target protospacer
    mismatches: int     # total mismatches vs the guide
    seed_mismatches: int
    risk_score: float   # 0.0 (low risk) → 1.0 (high risk)

    def __repr__(self):
        return (
            f"OffTargetHit(id={self.sequence_id}, pos={self.position}, "
            f"mm={self.mismatches}, seed_mm={self.seed_mismatches}, "
            f"risk={self.risk_score:.3f})"
        )


def _count_mismatches(spacer: str, candidate: str) -> tuple[int, int]:
    """
    Return (total_mismatches, seed_mismatches) between spacer and candidate.
    Both must be 20nt.
    """
    total = sum(1 for a, b in zip(spacer, candidate) if a != b)
    seed_mm = sum(
        1 for i in range(SEED_START, min(SEED_END, len(spacer)))
        if spacer[i] != candidate[i]
    )
    return total, seed_mm


def _mismatch_risk(spacer: str, candidate: str) -> float:
    """
    Weighted mismatch risk score (0.0 = identical, 1.0 = no shared seed).

    Seed mismatches are penalized more heavily than 5' mismatches because
    the seed region governs initial guide-target base pairing.
    """
    if len(spacer) != len(candidate):
        return 0.0
    weighted_penalty = 0.0
    for i in range(len(spacer)):
        if spacer[i] != candidate[i]:
            weight = MISMATCH_WEIGHTS.get(i, 0.5)
            weighted_penalty += weight

    max_penalty = sum(MISMATCH_WEIGHTS.get(i, 0.5) for i in range(len(spacer)))
    return round(1.0 - min(1.0, weighted_penalty / max_penalty), 4)


def find_offtargets(
    guide: GuideRNA,
    background_sequences: dict[str, str],
    max_mismatches: int = 3,
) -> list[OffTargetHit]:
    """
    Search background sequences for potential off-target cut sites.

    A site is flagged if it has an NGG PAM and the protospacer differs
    from the guide by at most max_mismatches.

    Args:
        guide: the guide RNA to check
        background_sequences: dict of {seq_id: sequence}
        max_mismatches: maximum mismatches to report (default 3)

    Returns:
        List of OffTargetHit sorted by risk score (highest first).
    """
    hits = []
    for seq_id, seq in background_sequences.items():
        pam_sites = find_pam_sites(seq)
        for site in pam_sites:
            if site.protospacer == guide.spacer:
                continue   # skip the on-target site itself
            total_mm, seed_mm = _count_mismatches(guide.spacer, site.protospacer)
            if total_mm <= max_mismatches:
                risk = _mismatch_risk(guide.spacer, site.protospacer)
                hits.append(OffTargetHit(
                    sequence_id=seq_id,
                    position=site.position,
                    strand=site.strand,
                    candidate=site.protospacer,
                    mismatches=total_mm,
                    seed_mismatches=seed_mm,
                    risk_score=risk,
                ))

    hits.sort(key=lambda h: h.risk_score, reverse=True)
    return hits


def offtarget_summary(hits: list[OffTargetHit]) -> dict:
    """Return aggregate off-target statistics for a guide."""
    if not hits:
        return {"total_hits": 0, "seed_hits": 0, "max_risk": 0.0, "mean_risk": 0.0}
    seed_hits = sum(1 for h in hits if h.seed_mismatches == 0)
    risks = [h.risk_score for h in hits]
    return {
        "total_hits": len(hits),
        "seed_hits": seed_hits,
        "max_risk": round(max(risks), 4),
        "mean_risk": round(sum(risks) / len(risks), 4),
    }
