#!/usr/bin/env python3
import argparse
import re
from Bio import SeqIO

TRANSCRIPT_SUFFIX = re.compile(r"\.t\d+$")

def parse_cds_fasta(fa_path):
    """
    解析OGAP cds fasta，提取 id, gene_name, cov
    返回 dict: gene_id(去掉.t1后缀) -> {"gene_name": str, "cov": float, "desc": str, "seq_record": SeqRecord}
    """
    gene_dict = {}
    for rec in SeqIO.parse(fa_path, "fasta"):
        desc = rec.description
        target_id = None
        gene_name = None
        cov = None
        for item in desc.split(";"):
            item = item.strip()
            if item.startswith("id="):
                target_id = item.split("=", 1)[1].strip()
            elif item.startswith("gene="):
                gene_name = item.split("=", 1)[1].strip()
            elif item.startswith("cov="):
                try:
                    cov = float(item.split("=", 1)[1].strip())
                except ValueError:
                    cov = None
        if target_id is None:
            continue
        # 裁剪 .t1/.t2 转录本后缀，对齐GFF gene ID
        gene_id = TRANSCRIPT_SUFFIX.sub("", target_id)
        gene_dict[gene_id] = {
            "gene_name": gene_name if gene_name is not None else gene_id,
            "cov": cov,
            "desc": desc,
            "seq_record": rec,
            "orig_transcript_id": target_id
        }
    return gene_dict

def parse_gff_genes(gff_path):
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
                    gid = kv.split("=",1)[1].strip()
                elif kv.startswith("Name="):
                    gname = kv.split("=",1)[1].strip()
            if gid:
                all_gene_ids.add(gid)
                gene_name_map[gid] = gname if gname else gid
    return gene_name_map, all_gene_ids

def filter_gff(in_gff, out_gff_path, bad_gene_ids):
    with open(in_gff,"r") as fin, open(out_gff_path,"w") as fout:
        for line in fin:
            keep = True
            if not line.startswith("#"):
                parts = line.split("\t")
                if len(parts)>=9:
                    ft = parts[2]
                    attr = parts[8]
                    gid = None
                    if ft == "gene":
                        for kv in attr.split(";"):
                            if kv.startswith("ID="):
                                gid = kv.split("=",1)[1].strip()
                                break
                    elif ft in ("CDS","exon"):
                        for kv in attr.split(";"):
                            if kv.startswith("Parent="):
                                gid = kv.split("=",1)[1].strip()
                                break
                    if gid is not None and gid in bad_gene_ids:
                        keep = False
            if keep:
                fout.write(line)

def main():
    parser = argparse.ArgumentParser(description="OGAP post-filter: filter gene by hmm coverage (cov=)")
    parser.add_argument("--gff", required=True, help="OGAP output gff3")
    parser.add_argument("--cds-fa", required=True, help="OGAP output cds fasta")
    parser.add_argument("--min-cov", type=float, default=80.0, help="minimum hmm coverage threshold, default 80.0")
    parser.add_argument("--mode", choices=["report","remove"], default="report")
    parser.add_argument("--out-gff", help="filtered gff output, required for remove mode")
    parser.add_argument("--out-cds", help="filtered cds fasta output, required for remove mode")
    parser.add_argument("--report", required=True, help="output tsv report")
    args = parser.parse_args()

    gene_dict = parse_cds_fasta(args.cds_fa)
    gff_gene_name_map, all_gene_ids = parse_gff_genes(args.gff)

    bad_gene_ids = set()
    report_rows = []

    for gid in all_gene_ids:
        gname = gff_gene_name_map.get(gid, gid)
        entry = gene_dict.get(gid)
        if entry is None:
            report_rows.append([gid, gname, "NA", "cds_sequence_missing"])
            continue
        cov = entry["cov"]
        if cov is None:
            status = "cov_missing"
        else:
            if cov < args.min_cov:
                status = "low_hmm_coverage"
                bad_gene_ids.add(gid)
            else:
                status = "ok"
        report_rows.append([gid, gname, str(cov if cov is not None else "NA"), status])

    with open(args.report,"w") as f:
        f.write("gene_id\tgene_name\tcov\tstatus\n")
        for row in report_rows:
            f.write("\t".join(row)+"\n")

    total = len(all_gene_ids)
    bad_cnt = len(bad_gene_ids)
    print(f"Total genes parsed: {total}")
    print(f"Genes with low_hmm_coverage: {bad_cnt}")

    if args.mode == "remove":
        if not args.out_gff or not args.out_cds:
            raise ValueError("--out-gff and --out-cds required under remove mode")
        filter_gff(args.gff, args.out_gff, bad_gene_ids)
        with open(args.out_cds,"w") as fout:
            for rec in SeqIO.parse(args.cds_fa,"fasta"):
                target_id = None
                for item in rec.description.split(";"):
                    if item.startswith("id="):
                        target_id = item.split("=",1)[1].strip()
                        break
                gid = None
                if target_id is not None:
                    gid = TRANSCRIPT_SUFFIX.sub("", target_id)
                if gid is not None and gid in bad_gene_ids:
                    continue
                fout.write(f">{rec.description}\n{rec.seq}\n")

if __name__ == "__main__":
    main()