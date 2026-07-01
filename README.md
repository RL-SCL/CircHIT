CircHIT – Quick Start & Reference
================================================================================

## 1. FOLDER STRUCTURE

After unzipping, you should see:

```text
CircHIT/
├── CircHIT.exe
├── _internal/
├── images/
├── 00_inputs/                 
│   ├── Raw_PFS_FD.csv                  (FD scores, required)
│   └── human_all_rna_v21_GRCh38.fa     (transcriptome FASTA, required)
├── 01_stage1_pfs_hits/                 (created automatically)
├── 02_stage2_scored/                   (created automatically)
└── 03_completed_logs/                  (created automatically)

```

> **Note:** The output folders appear automatically the first time you run an off‑target scan.

---

## 2. Sequence Inputs Guide

### Circular RNA Target

For a circular RNA, you only need two sequences (linear overhang and internal exons are no longer required):

* **BSJ Distal Sequence**: The exon *after* the back‑splice junction (BSJ).
* *Example:* For circMAP2K4(2‑4), this is exon 2.


* **BSJ Proximal Sequence** : The exon *before* the BSJ.
* *Example:* For circMAP2K4(2‑4), this is exon 4.



The program joins them **proximal → distal** to search the junction region.

### Fusion Gene Target

* **Fusion Distal Sequence:** The exon downstream of the fusion junction (e.g., the first exon of the *abl* transcript in the BCR‑ABL fusion).
* **Fusion Proximal Sequence:** The exon upstream of the fusion junction (e.g., the last exon of the *bcr* transcript).

The search logic is identical to the circular mode; only the column heading changes to **“Closest distance from Fusion”**.

### Gene (mRNA) Target

Simply paste the linear mRNA sequence (5′→3′, DNA letters). **Minimum 25 nt.**

---

## 3. Step‑by‑Step Workflow

### Circular RNA Target

1. Click **Circular RNA Target**.
2. Paste the DNA sequence of the **BSJ Distal exon** (minimum 25 nt). Click **Next**.
3. Paste the DNA sequence of the **BSJ Proximal exon** (minimum 25 nt). Click **Next**.
4. A list of candidate guideRNAs will appear in a new window.

#### Buttons in the Results Window

* **Export Results:** Saves the table as a CSV or Excel file.
* **Check Off‑Targets for Selected:** Select one row (click it), then click this button to scan the entire transcriptome for potential off‑target binding sites.
* **BLAST Selected Sequence:** Select a row and click to open NCBI BLASTN with the RC Target Region sequence pre‑filled. You can then choose a species, database, etc.

### Gene (mRNA) Target

1. Click **Gene (mRNA) Target**.
2. Paste the DNA sequence (5′→3′) of the mRNA (minimum 25 nt).
3. Click **Run Analysis**.
4. The results table appears with the same export, off‑target, and BLAST options.

### Fusion Gene Target

1. Click **Fusion Gene Target**.
2. Paste the DNA sequence of the **Fusion Distal exon** (minimum 25 nt). Click **Next**.
3. Paste the DNA sequence of the **Fusion Proximal exon** (minimum 25 nt). Click **Next**.
4. Results are displayed with the column “Closest distance from Fusion”.

### Batch Circular RNA Mode

1. Click **Batch Circular RNA Mode**.
2. Choose a CSV file with the following columns: `name`, `BSJ-Distal`, `BSJ-Proximal` (a template is provided in the download).
3. All rows are processed, and a combined table is displayed.

### Linear Batch Mode

1. Click **Linear Batch Mode**.
2. Choose a CSV file with columns: `name`, `sequence` (or at least `sequence`).
3. All sequences are analysed, and a combined table is displayed.

---

## 4. Interpreting Results

### A. Designer Output (All Modes)

| Column | Meaning |
| --- | --- |
| **PFS Motif** | The 5‑nt PFS motif identified. |
| **PFS Score** | Fold‑Depletion (FD) score. Higher = better. |
| **Closest distance from BSJ** (or **Fusion**) | Distance from the PFS to the junction. Negative = in the proximal exon, positive = in the distal exon. |
| **Target Region** | The 24‑nt window immediately upstream of the PFS. |
| **RC Target Region** | Reverse‑complement of the target region. **This is your candidate guide sequence (DNA).** |
| **guideRNA** | Full transcribed gRNA sequence: conserved tracrRNA handle + RC target region (RNA letters, U replaces T). |
| **DNA oligo for T7 IVT** | DNA oligo for ordering (target region + T7 promoter tail). |
| **GC % in sequence** | GC content of the RC target. |

### B. Off‑Target Scan Output

| Column | Meaning |
| --- | --- |
| **gene_name** | Gene or transcript identifier. |
| **rna_type** | RNA type (mRNA, tRNA, lncRNA, etc.). |
| **pfs_5mer** | PFS motif at the off‑target site. |
| **total_mismatches** | Total mismatches between guide and target. |
| **mismatch_intolerant** | Mismatches in the “seed” region (guide positions 1‑16). These pose a higher risk. |
| **mismatch_tolerant** | Mismatches in the distal region (17‑24). |
| **FD_score** | FD score of the off‑target PFS. |
| **protospacer_total_hits** | How many times the exact 24‑nt site appears in the transcriptome. |
| **protospacer_unique_records** | Unique transcripts containing that site. |
| **protospacer_unique_genes** | Unique genes containing that site. |
| **mismatch_positions_1based_all** | List of all mismatched positions (1‑based). |

> **Crucial Note:** Sites with few intolerant mismatches and low total mismatches present the highest off‑target risks.

---

## 5. First Off‑Target Scan

The first time you scan a guide, the program builds a large index of all PFS sites in the genome. This initial build may take **5–20 minutes**. All subsequent scans will be nearly instantaneous.

---

## 6. Troubleshooting

* **“FD CSV not found”:** Ensure `Raw_PFS_FD.csv` is located directly inside the `00_inputs` folder.
* **Scan seems stuck:** The background process is running. Please wait for the progress window to disappear.
* **“Missing columns in Stage‑1 TSV”:** Delete the `01_stage1_pfs_hits` folder and re‑run the scan to force the program to rebuild it.
