"""Formatted terminal output for guide RNA design results."""

from colorama import Fore, Style, init
init(autoreset=True)

DIVIDER = "─" * 64


def _header(title: str) -> str:
    return f"\n{Fore.CYAN}{Style.BRIGHT}{DIVIDER}\n  {title}\n{DIVIDER}{Style.RESET_ALL}"


def _label(key: str, value) -> str:
    return f"  {Fore.YELLOW}{key:<24}{Style.RESET_ALL}{value}"


def print_welcome():
    print(f"""
{Fore.CYAN}{Style.BRIGHT}
   ██████╗██████╗ ██╗███████╗██████╗ ██████╗
  ██╔════╝██╔══██╗██║██╔════╝██╔══██╗██╔══██╗
  ██║     ██████╔╝██║███████╗██████╔╝██████╔╝
  ██║     ██╔══██╗██║╚════██║██╔═══╝ ██╔══██╗
  ╚██████╗██║  ██║██║███████║██║     ██║  ██║
   ╚═════╝╚═╝  ╚═╝╚═╝╚══════╝╚═╝     ╚═╝  ╚═╝
{Style.RESET_ALL}
  {Fore.WHITE}CRISPR Guide RNA Designer — SpCas9 gRNA toolkit{Style.RESET_ALL}
""")


def print_pam_summary(seq_id: str, sites: list) -> None:
    print(_header(f"PAM Sites — {seq_id}"))
    fwd = [s for s in sites if s.strand == "+"]
    rev = [s for s in sites if s.strand == "-"]
    print(_label("Total PAM sites", len(sites)))
    print(_label("Forward strand (+)", len(fwd)))
    print(_label("Reverse strand (-)", len(rev)))
    print(f"\n  {'Pos':>5}  {'Str':^4}  {'PAM':^4}  Protospacer")
    print(f"  {'─'*5}  {'─'*4}  {'─'*4}  {'─'*20}")
    for s in sites[:15]:
        strand_color = Fore.GREEN if s.strand == "+" else Fore.MAGENTA
        print(
            f"  {s.position:>5}  "
            f"{strand_color}{s.strand:^4}{Style.RESET_ALL}  "
            f"{s.pam_seq:^4}  "
            f"{Fore.WHITE}{s.protospacer}{Style.RESET_ALL}"
        )
    if len(sites) > 15:
        print(f"  {Fore.YELLOW}... and {len(sites) - 15} more sites{Style.RESET_ALL}")


def print_ranked_guides(seq_id: str, ranked: list, top_n: int = 10) -> None:
    print(_header(f"Ranked Guide RNAs — {seq_id}"))
    if not ranked:
        print(f"  {Fore.RED}No viable guide RNAs found after filtering.{Style.RESET_ALL}")
        return
    print(_label("Candidates after GC filter", len(ranked)))
    print(f"\n  {'#':>3}  {'Score':^6}  {'GC%':^5}  {'Str':^4}  {'Pos':>5}  Spacer (5'→3')")
    print(f"  {'─'*3}  {'─'*6}  {'─'*5}  {'─'*4}  {'─'*5}  {'─'*20}")
    for rank, (guide, score) in enumerate(ranked[:top_n], 1):
        score_color = Fore.GREEN if score >= 0.65 else Fore.YELLOW if score >= 0.55 else Fore.RED
        strand_color = Fore.GREEN if guide.strand == "+" else Fore.MAGENTA
        flag_str = f"  {Fore.RED}⚠{Style.RESET_ALL}" if guide.flags else ""
        print(
            f"  {rank:>3}.  "
            f"{score_color}{score:.3f}{Style.RESET_ALL}  "
            f"{guide.gc_content*100:>4.0f}%  "
            f"{strand_color}{guide.strand:^4}{Style.RESET_ALL}  "
            f"{guide.position:>5}  "
            f"{Fore.WHITE}{guide.spacer}{Style.RESET_ALL}"
            f"{flag_str}"
        )
    if len(ranked) > top_n:
        print(f"  {Fore.YELLOW}... and {len(ranked) - top_n} more candidates{Style.RESET_ALL}")


def print_offtarget_report(guide_spacer: str, hits: list, summary: dict) -> None:
    print(_header(f"Off-Target Analysis — {guide_spacer[:12]}..."))
    print(_label("Off-target hits (≤3 mm)", summary["total_hits"]))
    print(_label("Seed region hits", summary["seed_hits"]))
    print(_label("Max risk score", f"{summary['max_risk']:.3f}"))
    print(_label("Mean risk score", f"{summary['mean_risk']:.3f}"))
    if hits:
        print(f"\n  {'Risk':^6}  {'MM':^4}  {'Seed MM':^8}  {'Str':^4}  {'Pos':>5}  Candidate")
        print(f"  {'─'*6}  {'─'*4}  {'─'*8}  {'─'*4}  {'─'*5}  {'─'*20}")
        for h in hits[:10]:
            risk_color = Fore.RED if h.risk_score >= 0.8 else Fore.YELLOW if h.risk_score >= 0.5 else Fore.GREEN
            print(
                f"  {risk_color}{h.risk_score:.3f}{Style.RESET_ALL}  "
                f"{h.mismatches:^4}  "
                f"{h.seed_mismatches:^8}  "
                f"{h.strand:^4}  "
                f"{h.position:>5}  "
                f"{Fore.WHITE}{h.candidate}{Style.RESET_ALL}"
            )
