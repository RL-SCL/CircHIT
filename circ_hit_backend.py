#!/usr/bin/env python3
"""
CircHIT Backend (updated)
-------------------------
Guide design (linear / circular) using FD scores and genome‑wide
off‑target scanner.

Column changes:
- "Term" removed.
- "Adjusted Term" → "PFS Motif"
- "gBlock" → "DNA oligo for IVT"

FD scores are now large integers (e.g., TAAGT=526) from the
new Raw_PFS_FD.csv file.
"""

import os
import csv
import gzip
import heapq
import itertools
import re
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

import pandas as pd

# ----------------------------------------------------------------------
# 1. Default motifs & gBlock tail
# ----------------------------------------------------------------------

DEFAULT_TERMS = [
    "TAGAC","TAGAA","TCAGT","TCAAC","TCAAT","TGATA","TGAGA","TGACA","TGAAC","TGATC","TGAGC",
    "TGACC","TGAAT","TGATT","TGAGT","TGACT","TGGAT","TTGAT","CAGAC","CAGAA","CCGAT","CCAGT",
    "CCAGG","CCAGC","CCACA","CCACT","CCACG","CCACC","CTAGC","CGGAT","CGATA","CGAGA","CGACA",
    "CGAAC","CGATC","CGAGC","CGACC","CGAAT","CGATT","CGAGT","CGACT","CACAT","CAGAC","CAGAA",
    "GAGAC","GAGAA","GGAAT","GGAGC","GGGAT","TAA","TAGAG","TATA","TACA","TCAGC","TCATG",
    "TCACA","TCGA","TGAAA","TGGAC","TTAA","TTATA","TTAGA","TTGAC","CAGT","CAGAG","CAA","CCGAC",
    "CCAGT","CCAGA","CTACA","CTAT","CTAGG","CTAA","CGGAC","CGAAA","CACAG","CATA","CAGT","CAGAG",
    "CAA","AAA","AAGA","AGAA","ACAA","ATAAT","GAA","GAGAG","GATA","GACA","GGAC","GACA","GGAC",
    "GGAAC","GGATC","GGAGA","GGGAC","GCAA","GCAGT","GCATG","GTAA"
]

DEFAULT_GBLOCK_TAIL = "ATCTACAATAGTAGAAATTTACCCTATAGTGAGTCGTATTACAG"

# ----------------------------------------------------------------------
# 2. FD score loading (unchanged – expects PFS,FD CSV)
# ----------------------------------------------------------------------

def load_fd_scores(csv_path: str) -> Dict[str, float]:
    """
    Reads a CSV with columns PFS,FD.
    Returns a dict mapping 5‑mer (upper case, T) → highest FD value.
    """
    fd_map: Dict[str, float] = {}
    if not csv_path or not os.path.exists(csv_path):
        print(f"[WARNING] FD CSV not found: {csv_path}. Using score 0 for all motifs.")
        return fd_map

    with open(csv_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        headers = [h.strip() for h in (reader.fieldnames or [])]
        if headers != ["PFS", "FD"]:
            raise ValueError(f"FD CSV must have columns 'PFS,FD'. Found: {headers}")
        for row in reader:
            pfs = (row.get("PFS") or "").strip().upper().replace("U", "T")
            if not pfs:
                continue
            try:
                fd_val = float(row.get("FD"))
            except (TypeError, ValueError):
                continue
            if pfs not in fd_map or fd_val > fd_map[pfs]:
                fd_map[pfs] = fd_val
    return fd_map


# ----------------------------------------------------------------------
# 3. Helper utilities
# ----------------------------------------------------------------------

def GC_calc(fa_string):
    _string = fa_string.upper()
    _G = _string.count('G')
    _C = _string.count('C')
    if len(_string) == 0:
        return 0.0
    return (_G + _C) / len(_string) * 100

def get_24nt_guide(rc_target: str) -> str:
    """
    Given a reverse‑complement target region (any length ≥ 24),
    return the 24‑nt guide by taking the PFS‑proximal portion (the last 24 nt).
    """
    if len(rc_target) < 24:
        raise ValueError("RC target must be at least 24 nt long")
    return rc_target[-24:]


# ----------------------------------------------------------------------
# 4. Guide designer (linear & circular) – with guideRNA, renamed oligo, column reorder
# ----------------------------------------------------------------------

def analyse_linear_sequence(sequence, terms=None, gblock_tail=None, fd_map=None):
    """
    Linear mode. 28‑nt target region.  Uses FD scores from fd_map.
    Output columns: PFS Motif, PFS Score, Target Region, RC Target Region,
                    guideRNA, DNA oligo for T7 IVT, GC % in sequence
    """
    sequence = sequence.upper()
    terms = terms if terms is not None else DEFAULT_TERMS
    gblock_tail = gblock_tail if gblock_tail is not None else DEFAULT_GBLOCK_TAIL
    if fd_map is None:
        fd_map = {}

    results = []
    for term in terms:
        for start_index in [i for i in range(len(sequence)) if sequence.startswith(term, i)]:
            upstream_start = max(0, start_index - 28)
            target_region = sequence[upstream_start:start_index]
            if len(target_region) != 28:
                continue

            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
            rc_target = "".join(complement.get(base, base) for base in reversed(target_region))
            gc_percent = round(GC_calc(rc_target), 1)

            adjusted_term = term
            if len(adjusted_term) < 5:
                pad_right = sequence[start_index + len(adjusted_term):start_index + 5]
                adjusted_term += pad_right

            pfs_score = fd_map.get(adjusted_term, 0.0)
            gblock_sequence = target_region + gblock_tail

            # Build guideRNA: conserved region + RC target converted to RNA (T→U)
            rna_rc_target = rc_target.replace('T', 'U')
            guide_rna = "UAAAUUUCUACUAUUGUAGAU" + rna_rc_target

            results.append({
                "PFS Motif": adjusted_term,
                "PFS Score": pfs_score,
                "Target Region": target_region,
                "RC Target Region": rc_target,
                "guideRNA": guide_rna,
                "DNA oligo for T7 IVT": gblock_sequence,
                "GC % in sequence": gc_percent,
            })

    df = pd.DataFrame(results)
    if df.empty:
        return df
    # Ensure column order
    desired_order = [
        "PFS Motif", "PFS Score", "Target Region", "RC Target Region",
        "guideRNA", "DNA oligo for T7 IVT", "GC % in sequence"
    ]
    df = df[desired_order]
    df = df.sort_values(by="PFS Score", ascending=False)
    df = df.drop_duplicates(subset=["Target Region"], keep="first")
    return df


def analyse_sequence(upstream, downstream, terms=None,
                     up_length=27, down_length=30, gblock_tail=None,
                     fd_map=None, min_dist=-7, max_dist=14, junction_label="BSJ"):   # default max_dist now 14
    """
    Circular mode. 24‑nt target region.  Uses FD scores.
    Returns (truncated_seq, truncated_upstream, truncated_downstream, bsj_site, df)

    Output columns:
        PFS Motif, PFS Score, Closest distance from BSJ,
        Target Region, RC Target Region, guideRNA,
        DNA oligo for T7 IVT, GC % in sequence
    """
    terms = terms if terms is not None else DEFAULT_TERMS
    gblock_tail = gblock_tail if gblock_tail is not None else DEFAULT_GBLOCK_TAIL
    if fd_map is None:
        fd_map = {}

    length_downstream = len(downstream)
    joint_sequence = downstream + upstream
    truncated_sequence = joint_sequence[length_downstream - down_length : length_downstream + up_length]
    truncated_upstream = upstream[:up_length]
    truncated_downstream = downstream[-17:]
    bsj_site = downstream[-1] + upstream[0]
    search_sequence = downstream[-4:] + upstream[:up_length]
    length_diff = len(joint_sequence) - len(upstream)

    results = []
    for term in terms:
        for start_index_trunc in [i for i in range(len(search_sequence)) if search_sequence.startswith(term, i)]:
            bsj_distance = start_index_trunc - 4
            if not (min_dist <= bsj_distance <= max_dist):
                continue

            start_index_long = start_index_trunc + length_diff
            target_start = max(0, start_index_long - 28)
            target_sequence = joint_sequence[target_start : start_index_long - 4]

            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
            rc_target_sequence = "".join(complement.get(base, base) for base in reversed(target_sequence))
            gc_percent = round(GC_calc(rc_target_sequence), 1)

            adjusted_sequence = term
            if len(adjusted_sequence) < 5:
                pad_right = search_sequence[start_index_trunc + len(adjusted_sequence):start_index_trunc + 5]
                adjusted_sequence += pad_right

            pfs_score = fd_map.get(adjusted_sequence, 0.0)
            gblock_sequence = target_sequence + gblock_tail

            # Build guideRNA: conserved region + RC target converted to RNA (T→U)
            rna_rc_target = rc_target_sequence.replace('T', 'U')
            guide_rna = "UAAAUUUCUACUAUUGUAGAU" + rna_rc_target

            results.append({
                "PFS Motif": adjusted_sequence,
                "PFS Score": pfs_score,
                f"Closest distance from {junction_label}": bsj_distance,
                "Target Region": target_sequence,
                "RC Target Region": rc_target_sequence,
                "guideRNA": guide_rna,
                "DNA oligo for T7 IVT": gblock_sequence,
                "GC % in sequence": gc_percent,
            })

    df = pd.DataFrame(results)
    if df.empty:
        return truncated_sequence, truncated_upstream, truncated_downstream, bsj_site, df

    desired_order = [
        "PFS Motif",
        "PFS Score",
        f"Closest distance from {junction_label}",
        "Target Region",
        "RC Target Region",
        "guideRNA",
        "DNA oligo for T7 IVT",
        "GC % in sequence",
    ]
    df = df[desired_order]
    df = df.sort_values(by="PFS Score", ascending=False)
    df = df.drop_duplicates(subset=[f"Closest distance from {junction_label}"], keep="first")
    return truncated_sequence, truncated_upstream, truncated_downstream, bsj_site, df
# ----------------------------------------------------------------------
# 5. Off‑target pipeline (Stage 1 + Stage 2)
# ----------------------------------------------------------------------

# IUPAC table & PFS set for Stage 1 scanning
RAW_PFS_TABLE = """
UAANN	CAANN	GAANN	AAANN
UAGAV	CAGAV	GAGAV	AAGAN
UAUAN	CAGUN	GAUAN
UACAN	CAUAN	GACAN
	CACAK

UCARY	CCASN	GCAUG	ACAAN
UCAUG	CCAGU	GCAGU
UCACA	CCGAY	GCAAN
UCGAN

UGANH	CGANH	GGACN	AGAAN
UGGAY	UGGAY	GGAAY
		GGAUC
		GGAGM
		GGGAY

UUAAN	CUAAN	GUAAN	AUAAU
UUAUA	CUAGS
UUAGA	CUAUN
UUGAY	CUACA
""".strip()

IUPAC_DNA = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "N": "ACGT", "R": "AG", "Y": "CT", "S": "GC",
    "W": "AT", "K": "GT", "M": "AC", "B": "CGT",
    "D": "AGT", "H": "ACT", "V": "ACG",
}

def dna(s: str) -> str:
    return (s or "").strip().upper().replace("U", "T")

def revcomp_dna(s: str) -> str:
    comp = str.maketrans({"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"})
    return dna(s).translate(comp)[::-1]

def open_fasta(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")

def normalize_seq(seq: str) -> str:
    seq = seq.upper().replace("U", "T")
    return re.sub(r"[^ACGTNRYSWKMBDHV]", "N", seq)

def fasta_reader(path: str) -> Iterable[Tuple[str, str]]:
    with open_fasta(path) as fh:
        header = None
        seq_chunks = []
        for line in fh:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if header is not None:
                    yield header, normalize_seq("".join(seq_chunks))
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, normalize_seq("".join(seq_chunks))

def parse_pfs_table(raw: str) -> List[str]:
    motifs = []
    for line in raw.splitlines():
        for cell in line.split():
            motifs.append(cell.upper().replace("U", "T"))
    return motifs

def expand_iupac(motif: str) -> List[str]:
    pools = [IUPAC_DNA[b] for b in motif]
    return ["".join(p) for p in itertools.product(*pools)]

def build_pfs_set() -> Set[str]:
    expanded = set()
    for motif in parse_pfs_table(RAW_PFS_TABLE):
        if len(motif) == 5:
            expanded.update(expand_iupac(motif))
    return expanded

def infer_rna_type(header: str) -> str:
    h = header.lower()
    if "trna" in h: return "tRNA"
    if "rrna" in h: return "rRNA"
    if "lncrna" in h: return "lncRNA"
    if "mirna" in h: return "miRNA"
    if "mrna" in h: return "mRNA"
    return "unknown"

def parse_ids(header: str) -> Tuple[str, str, str]:
    record_id = header.split()[0].strip()
    return record_id, record_id, record_id

def build_stage1(fasta_path: str, fd_csv: str, out_tsv: str, progress_callback=None):
    """
    Build the Stage‑1 hit TSV.  If progress_callback is provided, it will be
    called with (current_record, total_records) periodically.
    """
    pfs_set = build_pfs_set()
    print(f"[INFO] Expanded accepted PFS 5-mers: {len(pfs_set)}")
    fd_map = load_fd_scores(fd_csv)
    print(f"[INFO] FD scores loaded: {len(fd_map)}")

    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    with open(out_tsv, "w", encoding="utf-8", newline="") as out:
        out.write("record_id\tgene_id\tgene_name\trna_type\tpfs_5mer\tFD_score\tpfs_start_0\tupstream_24nt\n")
        record_count = 0
        # We don't know total records in advance, so we can't give a percent; just ping.
        for header, seq in fasta_reader(fasta_path):
            record_count += 1
            if progress_callback and record_count % 100 == 0:
                progress_callback(record_count, 0)  # total unknown
            rid, gene_id, gene_name = parse_ids(header)
            rna_type = infer_rna_type(header)
            for i in range(len(seq) - 4):
                pfs = seq[i:i+5]
                if pfs not in pfs_set:
                    continue
                if i < 24:
                    continue
                upstream = seq[i-24:i]
                fd = fd_map.get(pfs)
                fd_str = "" if fd is None else f"{fd:.6f}"
                out.write(f"{rid}\t{gene_id}\t{gene_name}\t{rna_type}\t{pfs}\t{fd_str}\t{i}\t{upstream}\n")
    print(f"[INFO] Stage 1 done: {record_count} records")


@dataclass
class MatchProfile:
    total_mismatches: int
    mismatch_intolerant: int
    mismatch_tolerant: int
    mismatch_positions_1based_all: List[int]
    mismatch_positions_1based_intolerant: List[int]
    mismatch_positions_1based_tolerant: List[int]

def profile_match(guide_24: str, aligned_target_24: str) -> MatchProfile:
    g = dna(guide_24)
    t = dna(aligned_target_24)
    if len(g) != 24 or len(t) != 24:
        raise ValueError("Guide and aligned target must be 24 nt.")
    pos_all, pos_intol, pos_tol = [], [], []
    for i in range(24):
        if g[i] != t[i]:
            pos_all.append(i+1)
            if i < 16:
                pos_intol.append(i+1)
            else:
                pos_tol.append(i+1)
    return MatchProfile(
        total_mismatches=len(pos_all),
        mismatch_intolerant=len(pos_intol),
        mismatch_tolerant=len(pos_tol),
        mismatch_positions_1based_all=pos_all,
        mismatch_positions_1based_intolerant=pos_intol,
        mismatch_positions_1based_tolerant=pos_tol,
    )

def build_sort_key(mp: MatchProfile, fd: Optional[float], aligned_target: str) -> Tuple:
    fd_sort = fd if fd is not None else -1.0
    return (-mp.mismatch_intolerant, -mp.total_mismatches, -mp.mismatch_tolerant, fd_sort, aligned_target)

def score_single_guide(in_tsv: str, guide: str, top_n: int = 1000, mm_cutoff: int = 5,
                       add_multiplicity: bool = True, add_mismatch_positions: bool = True) -> List[Dict]:
    """
    Scores a single 24‑nt guide against a Stage‑1 TSV.
    Returns a list of dicts (the top rows), sorted best‑first.
    """
    guide = dna(guide)
    if len(guide) != 24:
        raise ValueError("Guide must be 24 nt")
    heap = []
    counter = 0
    skipped_len = 0
    seen = 0
    required = {"record_id","gene_id","gene_name","rna_type","pfs_5mer","FD_score","pfs_start_0","upstream_24nt"}
    with open(in_tsv, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Missing columns in Stage‑1 TSV: {sorted(missing)}")
        for row in reader:
            seen += 1
            raw_target = dna(row.get("upstream_24nt", ""))
            if len(raw_target) != 24:
                skipped_len += 1
                continue
            aligned_target = revcomp_dna(raw_target)
            mp = profile_match(guide, aligned_target)
            if mp.total_mismatches >= mm_cutoff:
                continue
            fd_str = row.get("FD_score", "")
            fd = None
            if fd_str.strip():
                try:
                    fd = float(fd_str)
                except ValueError:
                    pass
            out = {
                "record_id": row["record_id"],
                "gene_id": row["gene_id"],
                "gene_name": row["gene_name"],
                "rna_type": row["rna_type"],
                "pfs_5mer": row["pfs_5mer"],
                "FD_score": "" if fd is None else f"{fd:.6f}",
                "pfs_start_0": row["pfs_start_0"],
                "upstream_24nt": raw_target,
                "aligned_target_24nt": aligned_target,
                "guide": guide,
                "mismatch_intolerant": str(mp.mismatch_intolerant),
                "mismatch_tolerant": str(mp.mismatch_tolerant),
                "total_mismatches": str(mp.total_mismatches),
            }
            if add_mismatch_positions:
                out["mismatch_positions_1based_all"] = ",".join(map(str, mp.mismatch_positions_1based_all))
                out["mismatch_positions_1based_intolerant"] = ",".join(map(str, mp.mismatch_positions_1based_intolerant))
                out["mismatch_positions_1based_tolerant"] = ",".join(map(str, mp.mismatch_positions_1based_tolerant))
            sort_key = build_sort_key(mp, fd, aligned_target)
            counter += 1
            item = (sort_key, counter, out)
            if len(heap) < top_n:
                heapq.heappush(heap, item)
            else:
                if item[0] > heap[0][0]:
                    heapq.heapreplace(heap, item)
    heap.sort(reverse=True)
    final_rows = [x[2] for x in heap]
    # Multiplicity second pass
    if add_multiplicity and final_rows:
        targets_needed = set(r["aligned_target_24nt"] for r in final_rows)
        total_hits = defaultdict(int)
        rec_sets = defaultdict(set)
        gene_sets = defaultdict(set)
        with open(in_tsv, newline="", encoding="utf-8") as fh2:
            reader2 = csv.DictReader(fh2, delimiter="\t")
            for row2 in reader2:
                raw = dna(row2.get("upstream_24nt",""))
                if len(raw)!=24:
                    continue
                aln = revcomp_dna(raw)
                if aln not in targets_needed:
                    continue
                total_hits[aln] += 1
                rid = row2.get("record_id","").strip()
                gid = row2.get("gene_id","").strip()
                if rid: rec_sets[aln].add(rid)
                if gid: gene_sets[aln].add(gid)
        for r in final_rows:
            at = r["aligned_target_24nt"]
            r["protospacer_total_hits"] = str(total_hits[at])
            r["protospacer_unique_records"] = str(len(rec_sets[at]))
            r["protospacer_unique_genes"] = str(len(gene_sets[at]))
    return final_rows


def run_off_target_for_guide(guide_24: str,
                             workdir: str,
                             top_n: int = 1000,
                             mm_cutoff: int = 5,
                             progress_callback=None) -> List[Dict]:
    """
    Convenience function that builds Stage‑1 if needed, then scores a single guide.
    progress_callback(status_message: str) can be called to report progress.
    """
    stage1_tsv = os.path.join(workdir, "01_stage1_pfs_hits", "stage_1_hits.tsv")
    fasta = os.path.join(workdir, "00_inputs", "human_all_rna_v21_GRCh38.fa")
    fd_csv = os.path.join(workdir, "00_inputs", "Raw_PFS_FD.csv")

    if not os.path.exists(stage1_tsv):
        if progress_callback:
            progress_callback("Building Stage‑1 database (this may take a while)...")
        def cb(n, total):
            if progress_callback:
                progress_callback(f"Stage‑1: processed {n} records...")
        build_stage1(fasta, fd_csv, stage1_tsv, progress_callback=cb)
    else:
        if progress_callback:
            progress_callback("Stage‑1 TSV found, skipping build.")

    if progress_callback:
        progress_callback(f"Scoring off‑targets for guide {guide_24}...")
    results = score_single_guide(stage1_tsv, guide_24, top_n=top_n, mm_cutoff=mm_cutoff)
    return results