# 🧬 DNA Sequence Analysis Toolkit

A lightweight, dependency-free Python toolkit for basic DNA sequence analysis.
Built as a bioinformatics portfolio project.

## Features
- FASTA file parsing (multi-record)
- Base composition and GC content
- Reverse complement
- Transcription (DNA → mRNA)
- Translation (DNA → protein) using the standard codon table
- ORF finder on both strands, all 3 reading frames

## Usage
```bash
python dna_toolkit.py sample.fasta
python dna_toolkit.py sample.fasta --min-orf 30
```

### Sample output
```
=== demo_gene_1 Synthetic test sequence ===
Length          : 83 bp
GC content      : 53.01%
ORFs (>= 20 bp) : 1 found
  strand + frame 1 1-24 (24 bp): MAIVMGR
```

## Run tests
```bash
pip install pytest
pytest
```

## Project structure
```
dna-toolkit/
├── dna_toolkit.py        # main module + CLI
├── test_dna_toolkit.py   # unit tests
├── sample.fasta          # demo input
└── requirements.txt
```

## Future ideas
- Protein molecular weight and isoelectric point
- Restriction site finder
- Pairwise sequence alignment (Needleman-Wunsch)
- Fetch sequences directly from NCBI

## Author
Vighnesh Gaikwad — MSc Bioinformatics
