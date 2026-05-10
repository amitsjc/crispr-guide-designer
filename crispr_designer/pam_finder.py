"""PAM site finder — locates all NGG sites for SpCas9 on both DNA strands."""

from dataclasses import dataclass

GUIDE_LENGTH = 20   # SpCas9 uses a 20nt spacer
PAM = "NGG"         # SpCas9 PAM on the non-template strand


@dataclass
class PAMSite:
    """A candidate PAM site and its associated protospacer region."""
    position: int       # 0-based position of the first PAM nucleotide
    strand: str         # "+" (forward) or "-" (reverse)
    pam_seq: str        # actual 3nt PAM (e.g. "AGG", "TGG")
    protospacer: str    # 20nt spacer immediately upstream of PAM (5'→3')

    def __repr__(self):
        return (
            f"PAMSite(pos={self.position}, strand={self.strand}, "
            f"PAM={self.pam_seq}, spacer={self.protospacer})"
        )


def _reverse_complement(seq: str) -> str:
    table = str.maketrans("ATGC", "TACG")
    return seq.upper().translate(table)[::-1]


def find_pam_sites(sequence: str) -> list[PAMSite]:
    """
    Find all SpCas9 PAM sites (NGG) in a DNA sequence on both strands.

    For the forward strand: protospacer is the 20nt immediately upstream of NGG.
    For the reverse strand: PAM is NCC on the forward sequence; the protospacer
    is read 5'→3' on the reverse complement strand.

    Returns sites sorted by position.
    """
    seq = sequence.upper()
    sites = []

    # Forward strand — look for xGG (any N followed by GG)
    for i in range(len(seq) - 2):
        if seq[i+1] == "G" and seq[i+2] == "G":
            pam_start = i
            spacer_start = pam_start - GUIDE_LENGTH
            if spacer_start < 0:
                continue
            protospacer = seq[spacer_start:pam_start]
            if len(protospacer) == GUIDE_LENGTH:
                sites.append(PAMSite(
                    position=pam_start,
                    strand="+",
                    pam_seq=seq[i:i+3],
                    protospacer=protospacer,
                ))

    # Reverse strand — NCC on forward = NGG on reverse complement
    for i in range(len(seq) - 2):
        if seq[i+1] == "C" and seq[i+2] == "C":
            # PAM (NCC) starts at i on forward; end of protospacer on fwd = i+3
            proto_end = i + 3 + GUIDE_LENGTH
            if proto_end > len(seq):
                continue
            fwd_proto = seq[i+3: proto_end]
            if len(fwd_proto) == GUIDE_LENGTH:
                protospacer = _reverse_complement(fwd_proto)
                sites.append(PAMSite(
                    position=i,
                    strand="-",
                    pam_seq=_reverse_complement(seq[i:i+3]),  # NGG on rev strand
                    protospacer=protospacer,
                ))

    sites.sort(key=lambda s: s.position)
    return sites
