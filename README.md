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
* The output folders appear the first time you run an off‑target scan.

## 2. WHAT THE SEQUENCE INPUTS MEAN (CIRCULAR RNA MODE)

For a circular RNA, you are asked for three sequence components:

* **Upstream sequence:** The part of the RNA *after* the back‑splice junction (BSJ). In a linear representation of the circle, this is the sequence that immediately follows the junction going 5'→3'.
* **Downstream sequence:** The part *before* the back‑splice junction. In the linear representation, this sits before the upstream sequence.
* **Linear overhang:** The region further upstream of the acceptor exon (i.e. extending 5' from the upstream sequence). It is used to assess similarity between the circular target and the linear pre‑mRNA, helping to avoid guides that also match the unspliced transcript.

**In Canonical Sequence:**
`5' – [linear overhang] – [upstream] – [internal exons] – [downstream] – 3'`

**In Linearised Circular Sequence:**
`5' – [downstream] – [BSJ] – [upstream] – [linear overhang] – 3'`

---

## 3. INTERPRETING RESULTS

### A. Designer Output (Circular / Linear)

| Column | Meaning |
| --- | --- |
| **Term** | The PFS motif that was found. |
| **Adjusted Term** | Motif padded to a full 5‑mer if needed. |
| **PFS Score** | Fold‑Depletion (FD) score. Higher = better. |
| **Within BSJ (circ only)** | True if the PFS is within 1 nt of the BSJ. |
| **Closest distance from BSJ** | Distance from the PFS to the BSJ. Negative values mean the motif starts in the downstream exon; positive = in the upstream. |
| **Target Region** | The 24‑nt window immediately upstream of the PFS. |
| **RC Target Region** | Reverse‑complement of the target region. **THIS IS YOUR CANDIDATE gRNA SEQUENCE.** |
| **GC % in sequence** | GC content of the RC target. |
| **Linear Overlap Similarity %** | How similar the target is to the linear overhang (circ only). Low values are safer. |
| **gBlock** | Target region + constant tail (for cloning). |

### B. Off‑Target Scan Output

| Column | Meaning |
| --- | --- |
| **gene_name** | Gene or transcript identifier. |
| **rna_type** | RNA type (mRNA, tRNA, lncRNA, etc.) |
| **pfs_5mer** | PFS motif at the off‑target site. |
| **total_mismatches** | Total mismatches between guide and target. |
| **mismatch_intolerant** | Mismatches in the "seed" region (guide positions 1‑16). These are more dangerous. |
| **mismatch_tolerant** | Mismatches in the distal region (17‑24). |
| **FD_score** | FD score of the off‑target PFS. |
| **protospacer_total_hits** | Times the exact 24‑nt site appears in the transcriptome. |
| **protospacer_unique_records** | Unique transcripts containing that site. |
| **protospacer_unique_genes** | Unique genes containing that site. |
| **mismatch_positions_1based_all** | List of all mismatched positions (1‑based). |

> **Crucial Note:** Sites with few intolerant mismatches and low total mismatches present the highest off‑target risks.

---

## 4. FIRST OFF‑TARGET SCAN

The first time you scan a guide, the program builds a large index of all PFS sites in the genome. This may take **5–20 minutes**. Subsequent scans will be nearly instant.

---

## 5. TROUBLESHOOTING

* **"FD CSV not found":** Make sure `Raw_PFS_FD.csv` is inside the `00_inputs` folder.
* **Scan seems stuck:** The process is running; please wait for the progress window to disappear.
* **"Missing columns in Stage‑1 TSV":** Delete the folder `01_stage1_pfs_hits` and re‑run the scan to rebuild it.



