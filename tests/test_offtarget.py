"""Tests for the off-target risk estimator."""

from crispr_designer.offtarget import (
    find_offtargets, offtarget_summary, _count_mismatches, _mismatch_risk,
    SEED_START
)
from crispr_designer.guide_designer import GuideRNA
from crispr_designer.pam_finder import PAMSite, GUIDE_LENGTH


def _make_guide(spacer: str) -> GuideRNA:
    site = PAMSite(position=20, strand="+", pam_seq="AGG", protospacer=spacer)
    return GuideRNA(spacer=spacer, pam_site=site, gc_content=0.50)


def _background(spacer: str, pam: str = "AGG") -> dict[str, str]:
    return {"bg": "A" * 5 + spacer + pam + "A" * 10}


# --- _count_mismatches ---

def test_identical_spacers_zero_mismatches():
    s = "ATGCATGCATGCATGCATGC"
    assert _count_mismatches(s, s) == (0, 0)

def test_one_mismatch_5prime():
    ref  = "ATGCATGCATGCATGCATGC"
    alt  = "TTGCATGCATGCATGCATGC"  # mismatch at pos 0 (outside seed)
    total, seed = _count_mismatches(ref, alt)
    assert total == 1
    assert seed == 0

def test_one_mismatch_in_seed():
    ref  = "ATGCATGCATGCATGCATGC"
    alt  = list(ref)
    alt[SEED_START] = "T" if ref[SEED_START] != "T" else "A"
    total, seed = _count_mismatches(ref, "".join(alt))
    assert total == 1
    assert seed == 1


# --- _mismatch_risk ---

def test_identical_max_risk():
    s = "ATGCATGCATGCATGCATGC"
    assert _mismatch_risk(s, s) == 1.0

def test_risk_decreases_with_mismatches():
    ref = "ATGCATGCATGCATGCATGC"
    one_mm = list(ref)
    one_mm[0] = "C"
    two_mm = list(ref)
    two_mm[0] = "C"; two_mm[1] = "G"
    r1 = _mismatch_risk(ref, "".join(one_mm))
    r2 = _mismatch_risk(ref, "".join(two_mm))
    assert r1 > r2

def test_seed_mismatch_penalizes_more():
    ref = "ATGCATGCATGCATGCATGC"
    fiveprime_mm = list(ref); fiveprime_mm[0] = "C"
    seed_mm = list(ref); seed_mm[SEED_START] = "C"
    r_five = _mismatch_risk(ref, "".join(fiveprime_mm))
    r_seed = _mismatch_risk(ref, "".join(seed_mm))
    assert r_seed < r_five   # seed mismatch → lower risk score (more penalized)


# --- find_offtargets ---

def test_no_hits_in_unrelated_sequence():
    guide = _make_guide("ATGCATGCATGCATGCATGC")
    background = {"bg": "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"}
    hits = find_offtargets(guide, background)
    assert hits == []

def test_exact_match_excluded_as_ontarget():
    spacer = "ATGCATGCATGCATGCATGC"
    guide = _make_guide(spacer)
    background = _background(spacer)
    # Exact match should be skipped (it's the on-target site)
    hits = find_offtargets(guide, background)
    assert not any(h.candidate == spacer for h in hits)

def test_near_match_within_threshold_detected():
    spacer = "ATGCATGCATGCATGCATGC"
    # Create a 2-mismatch variant
    near = list(spacer); near[0] = "C"; near[1] = "C"
    guide = _make_guide(spacer)
    background = _background("".join(near))
    hits = find_offtargets(guide, background, max_mismatches=3)
    assert len(hits) >= 1

def test_hits_sorted_by_risk_descending():
    spacer = "ATGCATGCATGCATGCATGC"
    guide = _make_guide(spacer)
    near1 = list(spacer); near1[0] = "C"
    near2 = list(spacer); near2[0] = "C"; near2[1] = "C"
    bg = "A" * 5 + "".join(near1) + "AGG" + "A" * 10
    bg += "A" * 5 + "".join(near2) + "TGG" + "A" * 10
    hits = find_offtargets(guide, {"bg": bg})
    if len(hits) >= 2:
        assert hits[0].risk_score >= hits[1].risk_score


# --- offtarget_summary ---

def test_summary_empty():
    s = offtarget_summary([])
    assert s["total_hits"] == 0
    assert s["max_risk"] == 0.0
