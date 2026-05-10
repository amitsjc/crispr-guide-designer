"""gRNA designer — filters PAM sites into viable guide RNAs."""

from dataclasses import dataclass, field
from crispr_designer.pam_finder import PAMSite, find_pam_sites

GC_MIN = 0.40   # below this, guide won't bind stably
GC_MAX = 0.70   # above this, risk of non-specific binding
HOMOPOLYMER_RUN = 4   # runs of 4+ identical bases reduce efficiency


@dataclass
class GuideRNA:
    """A filtered, viable CRISPR guide RNA candidate."""
    spacer: str          # 20nt targeting sequence (same as protospacer)
    pam_site: PAMSite
    gc_content: float    # fraction 0.0–1.0
    flags: list[str] = field(default_factory=list)  # any quality warnings

    @property
    def position(self) -> int:
        return self.pam_site.position

    @property
    def strand(self) -> str:
        return self.pam_site.strand

    def __repr__(self):
        gc_pct = f"{self.gc_content * 100:.0f}%"
        flag_str = f" ⚠ {', '.join(self.flags)}" if self.flags else ""
        return (
            f"GuideRNA(pos={self.position}, strand={self.strand}, "
            f"GC={gc_pct}, spacer={self.spacer}{flag_str})"
        )


def _gc_fraction(seq: str) -> float:
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def _homopolymer_run(seq: str) -> int:
    """Return the length of the longest run of a single nucleotide."""
    if not seq:
        return 0
    max_run = current_run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1
    return max_run


def _quality_flags(spacer: str) -> list[str]:
    """Return a list of quality warnings for a spacer sequence."""
    flags = []
    run = _homopolymer_run(spacer)
    if run >= HOMOPOLYMER_RUN:
        flags.append(f"homopolymer run of {run}")
    # Avoid poly-T (terminates RNA polymerase III transcription)
    if "TTTT" in spacer.upper():
        flags.append("poly-T stretch")
    return flags


def design_guides(
    sequence: str,
    gc_min: float = GC_MIN,
    gc_max: float = GC_MAX,
    exclude_flagged: bool = False,
) -> list[GuideRNA]:
    """
    Find all PAM sites and filter to viable guide RNAs.

    Filters applied:
      1. GC content must be within [gc_min, gc_max]
      2. Optionally exclude guides with quality flags (homopolymers, poly-T)

    Returns guides sorted by position.
    """
    pam_sites = find_pam_sites(sequence)
    guides = []

    for site in pam_sites:
        gc = _gc_fraction(site.protospacer)
        if not (gc_min <= gc <= gc_max):
            continue
        flags = _quality_flags(site.protospacer)
        if exclude_flagged and flags:
            continue
        guides.append(GuideRNA(
            spacer=site.protospacer,
            pam_site=site,
            gc_content=gc,
            flags=flags,
        ))

    return guides
