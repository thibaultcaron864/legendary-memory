#!/usr/bin/env python3
import sys
import os
import re
from statistics import mean

def parse_attributes(attr):
    """transform column 9 in dictionary"""
    d = {}
    for field in attr.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            d[k] = v
    return d

def avg_or_na(values):
    return f"{mean(values):.2f}" if values else "NA"

def main(gff_dir):
    print("\t".join([
        "Sample","Genes","mRNA","CDS","TRNA","RRNA",
        "Avg_gene_length(bp)","Avg_exons_per_gene","Avg_intron_length(bp)","Coding_percent",
        "CAZ_total","CAZ_GH","CAZ_AA","CAZ_CE","CAZ_CBM","CAZ_PL","CAZ_GT",
        "Secreted","Transmembrane","Protease",
        "PFAM_genes","InterPro_genes","EggNog_genes","COG_genes","GO_genes","antiSMASH_genes","BUSCO_genes",
        "Keywords_lignocellulose"
    ]))

    for fname in sorted(os.listdir(gff_dir)):
        if not fname.endswith((".gff", ".gff3")):
            continue
        sample = os.path.splitext(fname)[0]
        path = os.path.join(gff_dir, fname)

        genes, mrna, cds, trna, rrna = 0,0,0,0,0
        gene_lengths, cds_lengths = [], []
        exons_per_transcript, introns = [], []

        # sets for annotation
        pfam=set(); interpro=set(); eggnog=set(); cog=set(); go=set()
        antismash=set(); busco=set(); secreted=set(); tm=set(); protease=set()
        caz_total=set(); caz_GH=set(); caz_AA=set(); caz_CE=set(); caz_CBM=set(); caz_PL=set(); caz_GT=set()
        lignokey=set()

        # regex lignocellulose keywords
        ligno_pat = re.compile(r"(cellulase|cellobiohydrolase|endoglucanase|beta.?glucosidase|laccase|peroxidase|xylanase|xylosidase|arabinofuranosidase|hemicellulase|lignin)", re.I)

        transcripts = {}  # id -> list of (start,end)

        with open(path) as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.rstrip().split("\t")
                if len(parts) < 9: continue
                seqid, source, feature, start, end, score, strand, phase, attr = parts
                start, end = int(start), int(end)
                attrs = parse_attributes(attr)

                if feature == "gene":
                    genes += 1
                    gene_lengths.append(end - start + 1)
                elif feature == "mRNA":
                    mrna += 1
                elif feature == "CDS":
                    cds += 1
                    cds_lengths.append(end - start + 1)
                elif feature == "tRNA":
                    trna += 1
                elif feature == "rRNA":
                    rrna += 1

                # Exons: track for intron calc
                if feature == "exon":
                    parent = attrs.get("Parent","").split(",")[0]
                    parent = re.sub(r"[-.:].*","",parent)
                    transcripts.setdefault(parent, []).append((start,end))

                # Gene or transcript id for annotation
                gid = attrs.get("ID") or attrs.get("Parent","")
                if not gid: continue
                gid = re.sub(r"[-.:].*","",gid)

                # Column 9 for functional terms
                a = attr

                if "PFAM" in a.upper(): pfam.add(gid)
                if "INTERPRO" in a.upper(): interpro.add(gid)
                if "EGGNOG" in a.upper(): eggnog.add(gid)
                if re.search(r"COG[:=]",a,re.I): cog.add(gid)
                if re.search(r"GO[:=]",a,re.I): go.add(gid)
                if "ANTISMASH" in a.upper(): antismash.add(gid)
                if "BUSCO" in a.upper(): busco.add(gid)
                if re.search(r"SECRETED|SIGNALP",a,re.I): secreted.add(gid)
                if re.search(r"TRANSMEMBRANE|TMHELIX|TM_DOMAIN|TM[0-9]",a,re.I): tm.add(gid)
                if re.search(r"MEROPS|PROTEASE|PEPTIDASE",a,re.I): protease.add(gid)

                if "CAZy" in a:
                    caz_total.add(gid)
                    if re.search(r"GH[0-9]",a): caz_GH.add(gid)
                    if re.search(r"AA[0-9]",a): caz_AA.add(gid)
                    if re.search(r"CE[0-9]",a): caz_CE.add(gid)
                    if re.search(r"CBM[0-9]",a): caz_CBM.add(gid)
                    if re.search(r"PL[0-9]",a): caz_PL.add(gid)
                    if re.search(r"GT[0-9]",a): caz_GT.add(gid)

                if ligno_pat.search(a): lignokey.add(gid)

        # exons/introns stats
        for tx, exs in transcripts.items():
            exs.sort()
            exons_per_transcript.append(len(exs))
            for i in range(1,len(exs)):
                intron = exs[i][0] - exs[i-1][1] - 1
                if intron > 0: introns.append(intron)

        # averages
        avg_gene_len = avg_or_na(gene_lengths)
        avg_exons = avg_or_na(exons_per_transcript)
        avg_intron = avg_or_na(introns)
        coding_percent = f"{100*sum(cds_lengths)/sum(gene_lengths):.2f}" if gene_lengths else "NA"

        print("\t".join(map(str,[
            sample, genes, mrna, cds, trna, rrna,
            avg_gene_len, avg_exons, avg_intron, coding_percent,
            len(caz_total), len(caz_GH), len(caz_AA), len(caz_CE), len(caz_CBM), len(caz_PL), len(caz_GT),
            len(secreted), len(tm), len(protease),
            len(pfam), len(interpro), len(eggnog), len(cog), len(go), len(antismash), len(busco),
            len(lignokey)
        ])))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <gff_dir>", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1])

