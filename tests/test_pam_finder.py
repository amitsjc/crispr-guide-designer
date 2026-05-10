"""Tests for the PAM site finder."""

from crispr_designer.pam_finder import find_pam_sites, PAMSite, GUIDE_LENGTH


def test_finds_forward_pam():
    # 20nt spacer + AGG PAM
    spacer = "A" * GUIDE_LENGTH
    seq = spacer + "AGG"
    sites = find_pam_sites(seq)
    fwd = [s for s in sites if s.strand == "+"]
    assert len(fwd) == 1
    assert fwd[0].pam_seq == "AGG"
    assert fwd[0].protospacer == spacer


def test_finds_reverse_pam():
    # NCC on forward = NGG on reverse
    seq = "A" * 3 + "TCC" + "A" * GUIDE_LENGTH
    sites = find_pam_sites(seq)
    rev = [s for s in sites if s.strand == "-"]
    assert len(rev) >= 1


def test_no_pam_in_short_sequence():
    # Sequence too short to fit spacer + PAM
    seq = "ATGCATGCAGG"  # only 11 bp, spacer needs 20
    sites = find_pam_sites(seq)
    fwd = [s for s in sites if s.strand == "+"]
    assert fwd == []


def test_protospacer_length():
    spacer = "GCTAGCTAGCTAGCTAGCTA"  # 20nt
    seq = spacer + "CGG" + "AAAA"
    sites = find_pam_sites(seq)
    fwd = [s for s in sites if s.strand == "+"]
    assert all(len(s.protospacer) == GUIDE_LENGTH for s in fwd)


def test_multiple_pams_in_sequence():
    spacer = "A" * GUIDE_LENGTH
    seq = spacer + "AGG" + "N" * 5 + spacer + "TGG"
    sites = find_pam_sites(seq)
    fwd = [s for s in sites if s.strand == "+"]
    assert len(fwd) >= 2


def test_pam_seq_ends_with_gg():
    spacer = "G" * GUIDE_LENGTH
    seq = spacer + "CGG"
    sites = find_pam_sites(seq)
    fwd = [s for s in sites if s.strand == "+"]
    assert all(s.pam_seq.endswith("GG") for s in fwd)


def test_sorted_by_position():
    spacer = "A" * GUIDE_LENGTH
    seq = spacer + "AGG" + "NNNNNNNNNNNNNNNNNNNN" + "TGG"
    sites = find_pam_sites(seq)
    positions = [s.position for s in sites]
    assert positions == sorted(positions)


def test_empty_sequence():
    assert find_pam_sites("") == []
