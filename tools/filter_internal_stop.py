#!/usr/bin/env python3
import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq

TRANSCRIPT_SUFFIX = re.compile(r"\.t\d+$")

def load_cds_fasta(fa_path):
    """
    读取OGAP CDS fasta，解析id=xxx，裁剪.t1转录本后缀，返回 dict: gene_id -> Seq
    """
    cds_dict = {}
    for rec in SeqIO.parse(fa_path, "fasta"):
        desc = rec.description
        target_id = None
        for item in desc.split(";"):
            if item.startswith("id="):
                target_id = item.split("=", 1)[1].strip()
                break
        if target_id is None:
            target_id = rec.id
        gene_id = TRANSCRIPT_SUFFIX.sub("", target_id)
        cds_dict[gene_id] = rec.seq
    return cds_dict

def parse_gff_gene_map(gff_path):
    gene_name_map = {}
    all_gene_ids = set()
    with open(gff_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "gene":
                continue
            attr = parts[8]
            gid = None
            gname = None
            for kv in attr.split(";"):
                kv = kv.strip()
                if kv.startswith("ID="):
                    gid = kv.split("=", 1)[1].strip()
                elif kv.startswith("Name="):
                    gname = kv.split("=", 1)[1].strip()
            if gid:
                all_gene_ids.add(gid)
                gene_name_map[gid] = gname if gname else gid
    return gene_name_map, all_gene_ids

def filter_gff(in_gff, out_gff_path, bad_gene_ids):
    with open(in_gff, "r") as fin, open(out_gff_path, "w") as fout:
        for line in fin:
            keep = True
            if not line.startswith("#"):
                parts = line.split("\t")
                if len(parts) >= 9:
                    ft = parts[2]
                    attr = parts[8]
                    gid = None
                    if ft == "gene":
                        for kv in attr.split(";"):
                            if kv.startswith("ID="):
                                gid = kv.split("=", 1)[1].strip()
                                break
                    elif ft in ("CDS", "exon"):
                        for kv in attr.split(";"):
                            if kv.startswith("Parent="):
                                gid = kv.split("=", 1)[1].strip()
                                break
                    if gid is not None and gid in bad_gene_ids:
                        keep = False
            if keep:
                fout.write(line)

def main():
    parser = argparse.ArgumentParser(description="OGAP post-filter: detect internal stop codon in CDS")
    parser.add_argument("--gff", required=True, help="OGAP output gff3")
    parser.add_argument("--cds-fa", required=True, help="OGAP output cds fasta")
    parser.add_argument("--code", type=int, default=1, help="genetic code table number, default 1 standard")
    parser.add_argument("--mode", choices=["report", "remove"], default="report")
    parser.add_argument("--out-gff", help="filtered gff output, required for remove mode")
    parser.add_argument("--out-cds", help="filtered cds fasta output, required for remove mode")
    parser.add_argument("--report", required=True, help="output tsv report file")
    args = parser.parse_args()

    cds_dict = load_cds_fasta(args.cds_fa)
    gene_name_map, all_gene_ids = parse_gff_gene_map(args.gff)

    bad_gene_ids = set()
    report_rows = []

    for gid in all_gene_ids:
        gname = gene_name_map.get(gid, gid)
        if gid not in cds_dict:
            report_rows.append([gid, gname, "cds_sequence_missing", ""])
            continue
        dna_seq = cds_dict[gid]
        prot = dna_seq.translate(table=args.code, to_stop=False)
        prot_str = str(prot)
        internal_stop_pos = None
        # 只检查除末尾外的内部终止
        for idx, aa in enumerate(prot_str[:-1]):
            if aa == "*":
                internal_stop_pos = idx + 1
                break
        if internal_stop_pos is not None:
            status = "internal_stop"
            bad_gene_ids.add(gid)
        else:
            status = "ok"
        report_rows.append([gid, gname, status, str(internal_stop_pos if internal_stop_pos else "")])

    with open(args.report, "w") as f:
        f.write("gene_id\tgene_name\tstatus\tstop_aa_pos\n")
        for row in report_rows:
            f.write("\t".join(row) + "\n")

    total = len(all_gene_ids)
    bad_cnt = len(bad_gene_ids)
    print(f"Total genes parsed: {total}")
    print(f"Genes with internal stop codon: {bad_cnt}")

    if args.mode == "remove":
        if not args.out_gff or not args.out_cds:
            raise ValueError("--out-gff and --out-cds required under remove mode")
        filter_gff(args.gff, args.out_gff, bad_gene_ids)
        with open(args.out_cds, "w") as fout:
            for rec in SeqIO.parse(args.cds_fa, "fasta"):
                target_id = None
                for item in rec.description.split(";"):
                    if item.startswith("id="):
                        target_id = item.split("=", 1)[1].strip()
                        break
                gid = None
                if target_id is not None:
                    gid = TRANSCRIPT_SUFFIX.sub("", target_id)
                if gid is not None and gid in bad_gene_ids:
                    continue
                fout.write(f">{rec.description}\n{rec.seq}\n")

if __name__ == "__main__":
    main()