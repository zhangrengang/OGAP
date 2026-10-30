# OGAP — Organelle Genome Annotation Pipeline

> Note: Adapted for major Python3 upgrade.

## 1. Project Introduction

- **Tool purpose**: Automatic annotation pipeline for chloroplast (pt) / mitochondrial (mt) organelle genomes.
- **Applicable taxa**: Mainly plants; good performance for animal protein‑coding gene annotation.
- **Main outputs**:
  `cov_report_Ara.tsv`, `filtered.gff3`, `sequence.fasta.fsa`, `sequence.fasta.gff3`, `sequence.fasta.rna.fasta`, `filtered.cds.fasta`, `sequence.fasta.cds.fasta`, `sequence.fasta.gb`, `sequence.fasta.pep.fasta`, `stop_report.tsv`
  > Note: Some files are only generated under specific workflows.

- **Current version**: Python3 refactored release (legacy Python2 version archived in branch `old_python2`, no longer maintained).

---

## 2. Major Release Notes [Core section for this upgrade]

### 2.1 Underlying migration: Python2 → Python3.8+

- Python2 is no longer supported.
- Fixed compatibility issues for syntax, encoding, iterators and I/O.
- Compatibility adapted for external software: AUGUSTUS / HMMER / tRNAscan‑SE / BLAT / Exonerate.

### 2.2 Pipeline structural refactoring

- New standalone database‑building script `ogap_build_db.py` for building custom databases.
- Separate management for pt and mt databases; independent `hmms`, `rna.info`, `cds.info`.
- Supports automatic generation of AUGUSTUS training profile (prfl). Models produced by custom training can be used for annotation. The prfl‑generation step is affected by system `libstdc++`. Compatibility anomalies occur after upgrading `libstdc++` from 6.027 to 6.035; nevertheless, trained models remain functional for annotation. Users need to manually validate output prfl files.

### 2.3 New features

- **rRNA alias whitelist normalization module**: Standardize rRNA gene names and resolve duplicate‑counting artifacts caused by inconsistent naming in legacy reference databases.

### 2.4 Deprecations & adjustments

- Archive all original Python2 scripts in branch `old_python2`.
- Adjust partial default parameters for Exonerate and remove unstable logic.

---

## 3. Environment & Dependencies

### 3.1 Python Environment

> Conda environment is recommended to avoid conflicts with system‑level Python libraries.

- **Python**: 3.12.13; minimum compatible version 3.8, **Python2 is unsupported**.
- **Biopython**: 1.87; minimum version 1.85 if OrthoFinder extended functions are enabled.
- **networkx**: 3.6.1.

### 3.2 External Bioinformatics Software

**✅ Required for main annotation workflow**

- **blat**: 36x series recommended; mainstream releases are generally compatible.
- **exonerate**: 2.4.0; backward‑compatible, baseline support for 2.2.x.
- **hmmsearch (HMMER)**: 3.4; backward‑compatible with 3.1, 3.2, 3.3 series.
- **tRNAscan‑SE**: Version 2.0 or higher is recommended.
- **AUGUSTUS**: Baseline version 3.1; tested up to 3.5, main annotation workflow runs correctly.
  > ⚠️ AUGUSTUS note: Custom‑training workflow remains available, trained models can be used for annotation. prfl generation is affected by system `libstdc++`. After upgrading `libstdc++` from 6.027 to 6.035, prfl format anomalies may occur. Trained models still work for annotation; manually validate prfl output when using custom training.

**✅ Invoked only by auxiliary scripts (not used in main annotation)**

- **hmmbuild**: Used exclusively in `ogap_build_db.py` for database construction; part of the HMMER suite, keep version consistent with hmmsearch.

**⚠️ Extended tools (not distributed with OGAP, download manually if needed)**

- **orthofinder‑3.1.5**: Required only for homology‑related extended functions during database building. Not required for standard pt/mt annotation. Version **3.1.5** is tested and functional. OGAP does not ship this binary; obtain it manually if needed.

> Remaining tools: aragorn, asn2gb; no invocation in pipeline code, not required.

### 3.3 Non‑core but pipeline‑supporting functions & dependencies (install on demand)

> The functional modules below are actually invoked by OGAP and contribute to full‑pipeline operation. They are **not part of the main annotation workflow**, are not bundled with OGAP, and have no enforced version constraints. Install corresponding tools manually only when activating these features.

| Functional module | Dependent tools | Trigger |
|---|---|---|
| Automatic taxonomic mapping | ete3 | `-taxon` omitted; taxon inferred from organism string |
| tRNA secondary‑structure plotting | RNAplot (ViennaRNA), ps2pdf | `-trn_struct` |
| Genome map drawing | OGDraw, latex | `-draw_map` |
| Comparative genome map (overlay raw GenBank record) | OGDraw, latex | `-compare_map` |
| Repeat‑region annotation | vmatch, trf | `-repeat` |
| Multi‑annotation comparison (phylogenetics) | mafft, trimal, iqtree | `Comparative.py` |
| Multi‑annotation comparison (KaKs calculation) | ParaAT, KaKs_Calculator | `Comparative.py` |

---

## 4. Database Preparation (New key feature)

### 4.1 Database design logic

- Databases are partitioned by organelle type (pt/mt) and taxonomic group.
- Database file description: `hmms`, `rna.info`, `cds.info`, reference sequences.

### 4.2 Example command for database construction

```bash
python3 ./OGAP_master/ogap_build_db.py \
-organ pt \
-gbfiles ./*.gb \
-custom name \
--dbroot-dir ./bd_made/ \
--orthofinder-exe ./orthofinder‑3.1.5/orthofinder.py \
--overwrite
```

> Note: `--orthofinder‑exe` is an argument for homology‑related extended functions. It can be removed for ordinary database construction.

---

## 5. Quick Start

> Short parameters `-pt` / `-mt` are aliases for `--ogtype pt` / `--ogtype mt`.

### 5.1 Minimal running examples

**Chloroplast (pt) annotation example**

```bash
python3 ./OGAP_master/OGAP.py \
./sequence.fasta \
-pt \
-o ./test1/ \
-sp Arabidopsis_thaliana \
-taxon Brassicaceae
```

**Mitochondrial (mt) annotation example**

```bash
python3 ./OGAP_master/OGAP.py \
./sequence.fasta \
-mt \
-o ./test1/ \
-sp Arabidopsis_thaliana \
-taxon Magnoliopsida
```

**Locusta migratoria mitochondrial example (real‑world test command)**

```bash
python3 ./OGAP_master/OGAP.py \
./sequence.fasta \
-mt \
-min_rrn_score=0.6 \
-trn-opts='-O -X 5' \
-o ./test/ \
-sp Locusta_migratoria \
-taxon Insecta
```

### 5.2 Core parameter description

| Parameter | Description |
|---|---|
| `-pt` / `-mt` | Specify organelle type; alias for `--ogtype pt/mt` |
| `-sp` | Species name, matches AUGUSTUS species model |
| `-taxon` | Taxonomic group, used for internal pipeline filtering and scoring |
| `--no‑rrna‑alias` | Disable rRNA alias normalization, revert to original annotation baseline |
| `-min_rrn_score=0.6` | HMM score threshold for rRNA; 0.6 commonly‑used for mitochondria, adjust per taxon |
| `-trn-opts='‑O ‑X 5'` | Arguments passed to tRNAscan‑SE; commonly‑used settings for mitochondrial mode |
| `-o` | Output directory path |

---

## 6. Output File Description

### 6.1 List of major output files

- **GFF annotation file**
- **Extracted gene FASTA sequences**
- **Runtime summary logs** (focus on `summary by source` / `summary by gene type`)
- **Hint files** (AUGUSTUS‑compatible annotation hints)

### 6.2 Log interpretation notes

> Distinguish two concepts:
>
> - **Biological true multi‑copy**: e.g. duplicated rRNA genes in chloroplast IR regions, retained by the pipeline.
> - **Alias‑caused duplicate counting**: one gene with multiple names. This artifact is resolved by the rRNA alias normalization module.

---

## 7. Detailed Feature Description (for advanced users)

### 7.1 rRNA alias normalization module

- **Background**: Legacy reference databases contain mixed rRNA aliases (e.g. `rrn16` / `rrn16S`). Identical rRNA sequences with different names are counted as separate genes.
- **Mechanism**: Whitelist‑based mapping; only predefined aliases are converted.
- **Mapping example**:

```python
rrn_alias = {
    # chloroplast
    "S5": "rrn5",
    "rrn5S": "rrn5",
    "5S": "rrn5",
    "4.5S": "rrn4.5",
    "16S": "rrn16",
    "rrn16S": "rrn16",
    "23S": "rrn23",
    "rrn23S": "rrn23",
    # insect mitochondria
    "rrn12S": "rrnS",
    "rrn16S": "rrnL"
}
```

- **Runtime log marker**: Log prints `rRNA rename: raw=xxx, mapped=xxx` for alias‑conversion traceability.
- **Boundary note**: Only modifies gene name field. **Does not change prediction sensitivity, coordinates or sequences**.
- **Limitation**: Whitelist‑driven. Manually extend the dictionary when encountering unrecorded new rRNA aliases.

---

## 8. Known Issues & Limitations

- **tRNAscan‑SE**: Limited sensitivity for mitochondrial tRNA prediction in certain taxa; false negatives may occur.
- **rRNA aliases**: Whitelist‑driven. Unrecorded new rRNA aliases require manual dictionary extension.
- **AUGUSTUS prfl generation**: Affected by system `libstdc++`. Format anomalies appear after upgrading `libstdc++` from 6.027 to 6.035. **Main gene‑prediction workflow is unaffected. Custom‑trained models remain usable for annotation; users must manually validate prfl outputs**.

---

## 9. Citation, Contact & License

> To‑be‑filled later.
