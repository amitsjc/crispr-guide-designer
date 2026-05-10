"""CLI — two subcommands: design and pamsites."""

import argparse
import sys
from crispr_designer.pam_finder import find_pam_sites
from crispr_designer.guide_designer import design_guides
from crispr_designer.scorer import rank_guides
from crispr_designer.offtarget import find_offtargets, offtarget_summary
from crispr_designer import report


def _load_fasta(path: str) -> dict[str, str]:
    """Return {seq_id: sequence} for all records in a FASTA file."""
    sequences = {}
    current_id = None
    parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(parts).upper()
                    parts = []
                current_id = line[1:].split()[0]
            else:
                parts.append(line)
    if current_id:
        sequences[current_id] = "".join(parts).upper()
    return sequences


def cmd_design(args):
    report.print_welcome()
    targets = _load_fasta(args.fasta)
    background = _load_fasta(args.offtarget) if args.offtarget else {}

    for seq_id, seq in targets.items():
        print(f"\n  Target: {Fore_cyan(seq_id)}  ({len(seq)} bp)")
        guides = design_guides(seq, gc_min=args.gc_min, gc_max=args.gc_max,
                               exclude_flagged=args.strict)
        ranked = rank_guides(guides)
        report.print_ranked_guides(seq_id, ranked, top_n=args.top)

        if background and ranked:
            top_guide, _ = ranked[0]
            hits = find_offtargets(top_guide, background, max_mismatches=args.max_mm)
            summary = offtarget_summary(hits)
            report.print_offtarget_report(top_guide.spacer, hits, summary)
    print()


def cmd_pamsites(args):
    targets = _load_fasta(args.fasta)
    for seq_id, seq in targets.items():
        sites = find_pam_sites(seq)
        report.print_pam_summary(seq_id, sites)
    print()


def Fore_cyan(text):
    from colorama import Fore, Style
    return f"{Fore.CYAN}{text}{Style.RESET_ALL}"


def main():
    parser = argparse.ArgumentParser(
        prog="crispr_designer",
        description="CRISPR Guide RNA Designer — SpCas9 gRNA toolkit",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # design
    p_design = sub.add_parser("design", help="Design and rank guide RNAs for a target")
    p_design.add_argument("fasta", help="Target sequence FASTA file")
    p_design.add_argument("--offtarget", metavar="FASTA",
                          help="Background genome FASTA for off-target analysis")
    p_design.add_argument("--gc-min", type=float, default=0.40, metavar="F",
                          help="Minimum GC content (default: 0.40)")
    p_design.add_argument("--gc-max", type=float, default=0.70, metavar="F",
                          help="Maximum GC content (default: 0.70)")
    p_design.add_argument("--top", type=int, default=10, metavar="N",
                          help="Show top N guides (default: 10)")
    p_design.add_argument("--max-mm", type=int, default=3, metavar="N",
                          help="Max mismatches for off-target search (default: 3)")
    p_design.add_argument("--strict", action="store_true",
                          help="Exclude guides with quality flags")
    p_design.set_defaults(func=cmd_design)

    # pamsites
    p_pam = sub.add_parser("pamsites", help="List all PAM sites in a sequence")
    p_pam.add_argument("fasta", help="Input FASTA file")
    p_pam.set_defaults(func=cmd_pamsites)

    args = parser.parse_args()
    try:
        args.func(args)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
