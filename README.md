# CRISPR Guide RNA Designer

A Python command-line tool for designing and scoring CRISPR-Cas9 guide RNAs. Given a target DNA sequence, this tool identifies all candidate cut sites, scores guide RNAs for on-target activity, and estimates off-target risk — the core computational pipeline used in CRISPR research.

## Background

CRISPR-Cas9 gene editing requires a **guide RNA (gRNA)** — a 20-nucleotide sequence that directs the Cas9 protein to a specific location in the genome. The Cas9 protein then cuts both DNA strands, enabling gene knockout, correction, or insertion.

For Cas9 to cut, a **PAM sequence** (Protospacer Adjacent Motif) must be present immediately downstream of the target: `5'-NGG-3'` for the widely used *Streptococcus pyogenes* Cas9 (SpCas9). Choosing the best guide RNA requires balancing:

- **On-target activity** — will the guide efficiently direct cutting at the intended site?
- **Off-target risk** — could the guide accidentally cut elsewhere in the genome?
- **GC content** — optimal guides have 40–70% GC for stable binding
- **Sequence rules** — avoid homopolymer runs, prefer G at position 20

## Features

- **PAM Site Finder** — locate all NGG sites on both forward and reverse strands
- **gRNA Designer** — extract 20nt spacers, filter by GC content and sequence quality rules
- **On-Target Activity Scorer** — position-weight scoring based on published experimental data
- **Off-Target Risk Estimator** — seed region mismatch analysis (the primary safety concern in CRISPR therapy)
- **Ranked CLI Report** — run as a command-line tool, outputs a ranked guide table

## Requirements

- Python 3.8+
- See `requirements.txt`

## Installation

```bash
git clone https://github.com/amitsjc/crispr-guide-designer.git
cd crispr-guide-designer
pip install -r requirements.txt
```

## Usage

```bash
# Design guides for a target sequence
python -m crispr_designer design data/vegfa.fasta

# Score and rank guides with off-target analysis
python -m crispr_designer design data/vegfa.fasta --offtarget data/genome_sample.fasta

# Show all PAM sites (no filtering)
python -m crispr_designer pamsites data/vegfa.fasta
```

## Sample Targets

- `data/vegfa.fasta` — human VEGFA gene exon 1 (a classic CRISPR research target)
- `data/emx1.fasta` — human EMX1 gene (used in the landmark 2013 Cong et al. Science paper)
- `data/genome_sample.fasta` — short genome excerpt for off-target testing

## Project Structure

```
crispr-guide-designer/
├── crispr_designer/     # Core modules
├── data/                # Sample FASTA files
├── tests/               # Unit tests
├── requirements.txt
└── README.md
```

## Related Project

This tool is designed to work with sequences from the [DNA Sequence Analyzer](https://github.com/amitsjc/dna-sequence-analyzer) — FASTA files produced or analyzed there can be fed directly into this pipeline.
