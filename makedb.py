#!/usr/bin/env python3
"""
OGAP standalone database build utility
Build mt/pt classification sub‑database (hmms, prfl and related files)
Input: RefSeq genbank gbff/gbff.gz files.
"""
import argparse
import os
import sys
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from lib.Database import Database


def check_binary(name: str) -> bool:
    """Check whether external binary exists in PATH."""
    try:
        subprocess.run(
            [name, "-h"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False
        )
        return True
    except FileNotFoundError:
        return False

# ========== 在这里插入 OrthoFinder 相关代码 ==========
#OF_SCRIPT = "/media/nfs6/qxc8/tools1/ogap_tools1_forp3/orthofinder-3.1.5/orthofinder.py"

def get_of_subenv():
    env = os.environ.copy()
    path_list = env["PATH"].split(os.pathsep)
    new_path = [p for p in path_list if not p.endswith("/OGAP_master/bin")]
    env["PATH"] = os.pathsep.join(new_path)
    return env
# ======================================================

def main():
    parser = argparse.ArgumentParser(
        description="OGAP build mt/pt classification sub‑database"
    )
    parser.add_argument(
        "-gbfiles",
        nargs="+",
        required=True,
        help="RefSeq genbank gbff/gbff‑gz files"
    )
    parser.add_argument(
        "-organ",
        choices=["mt", "pt"],
        required=True,
        help="mt=mitochondrion; pt=plastid"
    )
    parser.add_argument(
        "--filter-taxon",
        type=str,
        default=None,
        help="Single taxon name (e.g. Solanaceae), filter input GenBank records by this taxon"
    )
    parser.add_argument(
        "--filter-taxon-list",
        type=str,
        default=None,
        help="Text file, one taxon per line for batch building, filter input GenBank records"
    )
    parser.add_argument(
        "-taxon",
        type=str,
        default=None,
        help="Custom database name; this name will be used later in OGAP with -taxon"
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Force overwrite existing target database directory"
    )
    parser.add_argument(
        "--dbroot-dir",
        type=str,
        default=None,
        help="Output root directory for built database. If omitted, use OGAP default: pt=./db/pt, mt=./db/mt (inside OGAP_master)."
    )
    parser.add_argument(
        "--ncpu",
        type=int,
        default=24,
        help="Number of cpu threads [default: %(default)s]"
    )
    parser.add_argument(
        "--orthofinder-exe",
        type=str,
        default=None,
        help="Path to OrthoFinder-3.x main script (OrthoFinder.py). If not provided, search PATH for orthofinder. OrthoFinder is NOT bundled with OGAP."
    )
    
    args = parser.parse_args()
    print(f"CMD PARSE: organ={args.organ}, taxon(dbname)={args.taxon}")

    # 互斥逻辑：-taxon（库名）不能搭配 --filter-taxon / --filter-taxon-list
    if args.taxon is not None:
        if args.filter_taxon or args.filter_taxon_list:
            print("ERROR: -taxon (database name) cannot be used together with --filter-taxon / --filter-taxon-list", file=sys.stderr)
            sys.exit(1)
        taxon_list = [None]  # 标记自定义库名模式，循环只跑一次
    else:
        # 筛选模式：必须提供 --filter-taxon / --filter-taxon-list，二者互斥
        if not args.filter_taxon and not args.filter_taxon_list:
            print("ERROR: must supply either --filter-taxon or --filter-taxon-list or -taxon", file=sys.stderr)
            sys.exit(1)
        if args.filter_taxon and args.filter_taxon_list:
            print("ERROR: cannot use --filter-taxon and --filter-taxon-list together", file=sys.stderr)
            sys.exit(1)
        taxon_list = []
        if args.filter_taxon:
            taxon_list.append(args.filter_taxon.strip())
        else:
            with open(args.filter_taxon_list, "r", encoding="utf‑8") as fh:
                for line in fh:
                    t = line.strip()
                    if t:
                        taxon_list.append(t)

    # validate input gbff files
    for gz_file in args.gbfiles:
        if not os.path.isfile(gz_file):
            print(f"ERROR: input file not found: {gz_file}", file=sys.stderr)
            sys.exit(1)

    # check required external tools
    for tool in ("hmmbuild", "augustus"):
        if not check_binary(tool):
            print(f"ERROR: required tool not found in PATH: {tool}", file=sys.stderr)
            sys.exit(1)

    # build each target
    for taxon in taxon_list:
        if args.taxon is not None:
            display_name = args.taxon
            print(f"\n===== Building {args.organ}‑{display_name} (custom mode, no taxon filter) =====")
        else:
            display_name = taxon
            print(f"\n===== Building {args.organ}‑{display_name} =====")

        # ===================== 核心：默认路径绑定 OGAP_master/db =====================
        if args.dbroot_dir is not None:
            final_dbroot = args.dbroot_dir
        else:
            final_dbroot = os.path.join(SCRIPT_DIR, "db")
            print(f"INFO: --dbroot-dir not provided, using default root: {final_dbroot}")
        # ============================================================================

        db_obj = Database(
            organ=args.organ,
            taxon=taxon,
            custom=args.taxon,
            dbrootdir=final_dbroot,
            gbfiles=args.gbfiles,
            ncpu=args.ncpu,
            cleanup=True,
            allow_single_taxon=(args.taxon is not None)
        )
        db_path = db_obj.dbdir

        if os.path.exists(db_path):
            if args.overwrite:
                print(f"Overwrite enabled, remove old directory: {db_path}")
                shutil.rmtree(db_path)
            else:
                print(f"SKIP: {db_path} already exists. Use --overwrite to rebuild.", file=sys.stderr)
                continue

        os.makedirs(db_path, exist_ok=True)
        try:
            # db_obj.makedb()
            db_obj.makedb(orthofinder_exe=args.orthofinder_exe)
            print(f"SUCCESS: database built at {db_path}")
        except Exception as exc:
            print(f"BUILD FAILED for {display_name}: {exc}", file=sys.stderr)
            if os.path.exists(db_path):
                shutil.rmtree(db_path)
                print(f"Clean incomplete database directory: {db_path}", file=sys.stderr)
            raise


if __name__ == "__main__":
    main()
