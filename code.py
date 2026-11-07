import argparse
from collections import Counter

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def read_fasta(path):
    """Return dict {header: sequence} from a FASTA file."""
    records, header, chunks = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    records[header] = "".join(chunks).upper()
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header:
        records[header] = "".join(chunks).upper()
    return records


def validate(seq):
    """Raise ValueError if seq has non-ACGT characters."""
    bad = set(seq) - set("ACGT")
    if bad:
        raise ValueError(f"Invalid DNA characters: {sorted(bad)}")


def gc_content(seq):
    """GC percentage of a DNA sequence."""
    return 100 * (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def reverse_complement(seq):
    return "".join(COMPLEMENT[b] for b in reversed(seq))


def transcribe(seq):
    """DNA (coding strand) -> mRNA."""
    return seq.replace("T", "U")


def translate(seq):
    """DNA -> protein, stops at the first stop codon."""
    protein = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODON_TABLE[seq[i:i + 3]]
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


def find_orfs(seq, min_len=30):
    """Find ORFs (ATG ... stop) on both strands. Returns list of dicts."""
    orfs = []
    for strand, s in (("+", seq), ("-", reverse_complement(seq))):
        for frame in range(3):
            i = frame
            while i < len(s) - 2:
                if s[i:i + 3] == "ATG":
                    for j in range(i, len(s) - 2, 3):
                        if CODON_TABLE[s[j:j + 3]] == "*":
                            if (j - i) >= min_len:
                                orfs.append({
                                    "strand": strand,
                                    "frame": frame + 1,
                                    "start": i + 1,
                                    "end": j + 3,
                                    "length_bp": j + 3 - i,
                                    "protein": translate(s[i:j + 3]),
                                })
                            i = j  # continue after this ORF
                            break
                i += 3
    return sorted(orfs, key=lambda o: -o["length_bp"])


def analyze(name, seq, min_orf):
    validate(seq)
    counts = Counter(seq)
    print(f"\n=== {name} ===")
    print(f"Length          : {len(seq)} bp")
    print(f"Base counts     : " + ", ".join(f"{b}={counts[b]}" for b in "ACGT"))
    print(f"GC content      : {gc_content(seq):.2f}%")
    print(f"Reverse comp.   : {reverse_complement(seq)[:60]}...")
    print(f"mRNA            : {transcribe(seq)[:60]}...")
    orfs = find_orfs(seq, min_orf)
    print(f"ORFs (>= {min_orf} bp) : {len(orfs)} found")
    for o in orfs[:3]:
        print(f"  strand {o['strand']} frame {o['frame']} "
              f"{o['start']}-{o['end']} ({o['length_bp']} bp): {o['protein'][:40]}")


def main():
    parser = argparse.ArgumentParser(description="DNA Sequence Analysis Toolkit")
    parser.add_argument("fasta", help="Path to FASTA file")
    parser.add_argument("--min-orf", type=int, default=30,
                        help="Minimum ORF length in bp (default: 30)")
    args = parser.parse_args()

    for name, seq in read_fasta(args.fasta).items():
        analyze(name, seq, args.min_orf)


if __name__ == "__main__":
    main()
