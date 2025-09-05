#!/usr/bin/env python3
import sys
import os
import gzip
from statistics import median

def read_fasta_lengths(fasta_path):
    """get a list of contigs/scaffolds lengths"""
    lengths = []
    seq_len = 0
    opener = gzip.open if fasta_path.endswith(".gz") else open
    with opener(fasta_path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if seq_len > 0:
                    lengths.append(seq_len)
                seq_len = 0
            else:
                seq_len += len(line.strip())
        if seq_len > 0:
            lengths.append(seq_len)
    return lengths

def compute_Nx(lengths, x):
    """return N values"""
    total = sum(lengths)
    threshold = total * (x / 100.0)
    s = 0
    for L in sorted(lengths, reverse=True):
        s += L
        if s >= threshold:
            return L
    return 0

def compute_Lx(lengths, x):
    """return L values ."""
    total = sum(lengths)
    threshold = total * (x / 100.0)
    s = 0
    for i, L in enumerate(sorted(lengths, reverse=True), 1):
        s += L
        if s >= threshold:
            return i
    return len(lengths)

def gc_content(fasta_path):
    """count %GC."""
    g = c = a = t = 0
    opener = gzip.open if fasta_path.endswith(".gz") else open
    with opener(fasta_path, "rt") as f:
        for line in f:
            if not line.startswith(">"):
                seq = line.strip().upper()
                g += seq.count("G")
                c += seq.count("C")
                a += seq.count("A")
                t += seq.count("T")
    return 100.0 * (g + c) / (a + t + g + c) if (a + t + g + c) > 0 else 0

def count_Ns(fasta_path):
    """Compte le nombre de N (gaps) dans le fasta."""
    n_count = 0
    opener = gzip.open if fasta_path.endswith(".gz") else open
    with opener(fasta_path, "rt") as f:
        for line in f:
            if not line.startswith(">"):
                n_count += line.upper().count("N")
    return n_count

def main(fasta_dir):
    print("\t".join([
        "Sample", "Scaffolds", "Genome_length(bp)", "N50", "N90", "L50", "L90",
        "Longest_contig", "Ns_count", "Ns_percent", "GC%"
    ]))

    for fname in sorted(os.listdir(fasta_dir)):
        if not fname.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")):
            continue
        sample = os.path.splitext(fname)[0].replace(".fa", "").replace(".fna", "").replace(".fasta", "")
        path = os.path.join(fasta_dir, fname)
        lengths = read_fasta_lengths(path)
        if not lengths:
            continue
        total_len = sum(lengths)
        N50 = compute_Nx(lengths, 50)
        N90 = compute_Nx(lengths, 90)
        L50 = compute_Lx(lengths, 50)
        L90 = compute_Lx(lengths, 90)
        longest = max(lengths)
        Ns = count_Ns(path)
        Ns_percent = 100.0 * Ns / total_len if total_len > 0 else 0
        GC = gc_content(path)

        print("\t".join(map(str, [
            sample, len(lengths), total_len, N50, N90, L50, L90, longest,
            Ns, f"{Ns_percent:.3f}", f"{GC:.2f}"
        ])))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <fasta_dir>", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1])

