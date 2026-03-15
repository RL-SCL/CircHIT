import pandas as pd

# Default list of motif terms
DEFAULT_TERMS = [
    "TAGAC", "TAGAA", "TCAGT", "TCAAC", "TCAAT", "TGATA", "TGAGA", "TGACA", "TGAAC", "TGATC", "TGAGC",
    "TGACC", "TGAAT", "TGATT", "TGAGT", "TGACT", "TGGAT", "TTGAT", "CAGAC", "CAGAA", "CCGAT", "CCAGT",
    "CCAGG", "CCAGC", "CCACA", "CCACT", "CCACG", "CCACC", "CTAGC", "CGGAT", "CGATA", "CGAGA", "CGACA",
    "CGAAC", "CGATC", "CGAGC", "CGACC", "CGAAT", "CGATT", "CGAGT", "CGACT", "CACAT", "CAGAC", "CAGAA",
    "GAGAC", "GAGAA", "GGAAT", "GGAGC", "GGGAT", "TAA", "TAGAG", "TATA", "TACA", "TCAGC", "TCATG",
    "TCACA", "TCGA", "TGAAA", "TGGAC", "TTAA", "TTATA", "TTAGA", "TTGAC", "CAGT", "CAGAG", "CAA", "CCGAC",
    "CCAGT", "CCAGA", "CTACA", "CTAT", "CTAGG", "CTAA", "CGGAC", "CGAAA", "CACAG", "CATA", "CAGT", "CAGAG",
    "CAA", "AAA", "AAGA", "AGAA", "ACAA", "ATAAT", "GAA", "GAGAG", "GATA", "GACA", "GGAC", "GACA", "GGAC",
    "GGAAC", "GGATC", "GGAGA", "GGGAC", "GCAA", "GCAGT", "GCATG", "GTAA"
]

# Global constant for gBlock tail
DEFAULT_GBLOCK_TAIL = "ATCTACAATAGTAGAAATTTACCCTATAGTGAGTCGTATTACAG"


# GC percent calculation
def GC_calc(fa_string):
    _string = fa_string.upper()
    _G = _string.count('G')
    _C = _string.count('C')
    if len(_string) == 0:
        return 0.0
    return (_G + _C) / len(_string) * 100


def compare_sequence_similarity(A: str, B: str) -> float:
    if len(B) > len(A):
        return 0.0  # Graceful fallback instead of raising ValueError
    if len(B) == 0:
        return 0.0

    truncated_A = A[-len(B):]
    matches = sum(1 for a, b in zip(truncated_A, B) if a == b)
    return (matches / len(B)) * 100


def get_published_pfs_score(kmer):
    """Returns the published PFS score, or 0.0 if not found."""
    scoring_table = {
        "TAAGT": 9.04, "TAAGC": 9.02, "TAAGA": 8.99, "CAACA": 8.86, "CAAAG": 8.73,
        "GAACA": 8.65, "CAGAT": 8.62, "CAAAT": 8.54, "TAAAT": 8.49, "GAAAC": 8.47,
        "TAAAA": 8.47, "GAAAT": 8.43, "TCAAT": 8.42, "CAAAC": 8.34, "CAATG": 8.33,
        "TAACG": 8.31, "CAGAA": 8.23, "TCAAA": 8.20, "GAAAA": 8.15, "CAACC": 8.09,
        "CATAG": 8.06, "TAATA": 8.05, "TAAAG": 8.04, "CAATC": 7.97, "TCAAC": 7.96,
        "TAATT": 7.94, "TAACC": 7.84, "GTAAT": 7.80, "TAGAT": 7.77, "TAATG": 7.76,
        "CAAAA": 7.75, "CATAT": 7.73, "TACAT": 7.73, "CTAAG": 7.70, "TTAAA": 7.70,
        "CGAAG": 7.63, "TATAA": 7.63, "CTAAT": 7.63, "TGAAA": 7.62, "TAGAC": 7.59,
        "GAATA": 7.54, "CACAG": 7.53, "AAAAA": 7.45, "CGAAC": 7.42, "CAACT": 7.39,
        "GAACC": 7.37, "TGAAT": 7.36, "GAAAG": 7.36, "TATAG": 7.31, "CATAC": 7.30,
        "GAATC": 7.29, "TAAAC": 7.28, "CCGAT": 7.28, "TTAAT": 7.27, "GAGAG": 7.27,
        "GATAT": 7.25, "CATAA": 7.24, "AAAAC": 7.24, "TTAAG": 7.23, "TGAAC": 7.18,
        "CAATA": 7.16, "GAGAC": 7.14, "GGAAT": 7.14, "CAGAG": 7.13, "GCAAA": 7.11,
        "TAATC": 7.05, "CAAGA": 7.03, "GCAAC": 7.00, "TACAA": 6.99, "GCAAT": 6.99,
        "TAGAG": 6.97, "CAACG": 6.96, "CTAAC": 6.95, "TATAT": 6.92, "CAATT": 6.88,
        "CAGAC": 6.88, "AAAAG": 6.81, "GAATT": 6.77, "TAACT": 6.74, "CACAC": 6.73,
        "GATAC": 6.73, "GAAGA": 6.71, "TAGTA": 6.70, "TAAGG": 6.68, "CACAT": 6.67,
        "GAATG": 6.65, "TACAC": 6.63, "CAAGT": 6.63, "CGAGA": 6.63, "TACAG": 6.61,
        "TCAAG": 6.59, "CGATC": 6.57, "GAAGC": 6.54, "GTAAA": 6.54, "TAGCC": 6.47,
        "GAGAT": 6.46, "GGAAC": 6.43, "GCAAG": 6.40, "TGACA": 6.40, "TTAAC": 6.39,
        "CACAA": 6.39, "GATAA": 6.38, "GCAGA": 6.37, "GGACC": 6.31, "TATAC": 6.31,
        "GTAAC": 6.31, "AGAAG": 6.31, "GTAAG": 6.31, "GGATC": 6.30, "TAACA": 6.28,
        "TGGAC": 6.27, "CGAAA": 6.26, "CAAGC": 6.25, "CAGTA": 6.22, "TGGAT": 6.22,
        "AGAAT": 6.19, "TAGAA": 6.19, "GACAC": 6.18, "CTAAA": 6.16, "CGATA": 6.14,
        "GATAG": 6.10, "TGAGA": 6.09, "TGATC": 6.07, "GAAGT": 6.06, "GACAG": 6.05,
        "CCAGA": 6.05, "CAAGG": 6.04, "CAGTC": 6.04, "CCATA": 6.03, "CGAAT": 6.02,
        "CGAGC": 5.97, "CAGTT": 5.91, "CGACT": 5.90, "AAAGA": 5.88, "GAACT": 5.84,
        "TGATA": 5.84, "CAGCT": 5.82, "TGACC": 5.82, "GGACA": 5.81, "CGGAC": 5.81,
        "TCGAC": 5.80, "CGACC": 5.80, "GACAT": 5.78, "TCGAT": 5.77, "TGACT": 5.76,
        "TAGTG": 5.74, "AAGAT": 5.74, "CGGAT": 5.73, "AAACA": 5.73, "TCGAA": 5.71,
        "CGATT": 5.71, "AGAAC": 5.71, "AAATA": 5.70, "TGAGC": 5.69, "TGATT": 5.64,
        "TGAGT": 5.62, "CCAAA": 5.61, "CCGAC": 5.57, "CTAGA": 5.57, "GGGAC": 5.56,
        "AAAAT": 5.55, "CTGAT": 5.54, "AGAAA": 5.53, "GGAGA": 5.50, "TCATG": 5.50,
        "CTACA": 5.48, "TCATA": 5.47, "GAACG": 5.47, "TTACA": 5.45, "CAGTG": 5.43,
        "TTAGA": 5.42, "CAGGC": 5.42, "TCACA": 5.39, "CAGCA": 5.38, "CTAGC": 5.38,
        "GGAGC": 5.37, "CAGGA": 5.36, "CTATT": 5.36, "CGACA": 5.35, "CCACA": 5.34,
        "CTGAC": 5.33, "CTATA": 5.33, "TGACG": 5.33, "TAGCA": 5.32, "CTACT": 5.32,
        "TAGGA": 5.31, "TTACT": 5.30, "CGAGG": 5.30, "TTATC": 5.29, "TCATT": 5.28,
        "CTATC": 5.27, "GACAA": 5.26, "TAGTT": 5.25, "TAGGG": 5.25, "TTGAT": 5.24,
        "CAGCC": 5.24, "TCAGA": 5.24, "CCATT": 5.23, "CCAAG": 5.23, "CTACC": 5.23,
        "TAGTC": 5.23, "GGACT": 5.22, "TTATG": 5.21, "GCAGT": 5.19, "TTAGT": 5.18,
        "TAGCT": 5.18, "CGAGT": 5.17, "TTGAC": 5.15, "TCATC": 5.12, "TTATT": 5.11,
        "CTACG": 5.10, "TCACT": 5.09, "AAACG": 5.08, "CCACT": 5.07, "CTAGG": 5.07,
        "ACAAT": 5.06, "AAAGC": 5.04, "GAGTA": 5.03, "TGTAT": 5.03, "AAATC": 5.01,
        "TCACG": 5.01, "CTAGT": 5.01, "CGTAT": 5.00, "AAATT": 4.99, "TGATG": 4.99,
        "CCAAC": 4.99, "GAGAA": 4.98, "GGACG": 4.96, "CCAAT": 4.95, "GAGCA": 4.94,
        "GCACA": 4.94, "GGAAA": 4.94, "TAGCG": 4.91, "CCAGT": 4.91, "GGATG": 4.89,
        "CCTAT": 4.89, "CCAGC": 4.89, "AAAGG": 4.88, "AAGAC": 4.88, "CCAGG": 4.87,
        "GCGAC": 4.86, "CGACG": 4.83, "AAGAG": 4.83, "GGGAT": 4.83, "CAGGG": 4.82,
        "CCATC": 4.82, "TTATA": 4.81, "TTACC": 4.81, "GAGCC": 4.80, "GGAGT": 4.80,
        "CTATG": 4.79, "AAACT": 4.78, "TGAGG": 4.77, "TTAGG": 4.77, "ACAAA": 4.75,
        "AAAGT": 4.74, "GCGAT": 4.73, "CCATG": 4.72, "GAGTT": 4.70, "CAGCG": 4.68,
        "GAGCT": 4.67, "AAATG": 4.66, "TCAGT": 4.64, "GTATA": 4.64, "ATAAT": 4.63,
        "GAGTG": 4.63, "AAGAA": 4.61, "GGATA": 4.60, "TAGGC": 4.60, "GTATT": 4.58,
        "TCACC": 4.57, "TTACG": 4.57, "TGTAC": 4.57, "GCATT": 4.56, "GAAGG": 4.56,
        "CCACC": 4.55, "TATCT": 4.54, "CTTAT": 4.54, "TGTAG": 4.53, "GCATG": 4.52,
        "TCTAT": 4.51, "TAGGT": 4.51, "TATCG": 4.49, "TTGAA": 4.49, "TCGAG": 4.47,
        "ACAAG": 4.47, "GAGTC": 4.46, "TCAGC": 4.45, "GCAGC": 4.44, "CCACG": 4.43,
        "GCATA": 4.43, "TATCA": 4.43, "TGAAG": 4.41, "GTAGA": 4.40, "GCAGG": 4.39,
        "TTGAG": 4.39, "GATCA": 4.38, "GTACA": 4.37, "TCAGG": 4.36, "TTAGC": 4.35,
        "GTATC": 4.35, "GCATC": 4.34, "GCGAG": 4.30, "AAACC": 4.29, "ACAAC": 4.27,
        "CCTAC": 4.26, "GTAGT": 4.26, "CATCA": 4.22, "CAGGT": 4.19, "GTGAT": 4.19,
        "TGTAA": 4.18, "GCACG": 4.16, "TATGA": 4.14, "GTACT": 4.14, "CGTAC": 4.13,
        "TCTAC": 4.12, "CATGT": 4.09, "CGATG": 4.08, "GGATT": 4.08, "GCACT": 4.07,
        "TGGAG": 4.03, "GAGCG": 4.02, "GTGAG": 3.97, "TATGT": 3.92, "AATAT": 3.92,
        "TCTAG": 3.89, "CTTAC": 3.89, "GCACC": 3.88, "CATCG": 3.87, "CCGAG": 3.83,
        "CTTAG": 3.82, "GTATG": 3.76, "GATCC": 3.74, "GTAGG": 3.73, "GAGGT": 3.73,
        "GGAGG": 3.72, "GGTAT": 3.72, "GATCG": 3.71, "TATGG": 3.70, "GTGAC": 3.68,
        "TATCC": 3.66, "GTACG": 3.64, "AATAA": 3.63, "TATTT": 3.61, "CATGG": 3.54,
        "TATTG": 3.53, "CATTT": 3.51, "CGTAG": 3.49, "CATTA": 3.49, "CGGAG": 3.49,
        "CCTAG": 3.49, "CATCC": 3.48, "GGTAC": 3.47, "CCCAT": 3.47, "AATAG": 3.47,
        "TACGA": 3.46, "CGTAA": 3.45, "CATTG": 3.43, "ATAAG": 3.42, "CATGA": 3.41,
        "TATTA": 3.34, "GAGGC": 3.33, "TCTAA": 3.32, "CATGC": 3.31, "TGCAC": 3.31,
        "TTTAT": 3.28, "TGCAT": 3.26, "CTCAT": 3.24, "GAGGA": 3.24, "CTTAA": 3.24,
        "GGTAG": 3.17, "TATGC": 3.17, "GAGGG": 3.12, "CATCT": 3.12, "GATCT": 3.01,
        "GTACC": 2.99, "AGATA": 2.94, "AACAT": 2.91, "CTGAG": 2.91, "CATTC": 2.87,
        "GTAGC": 2.84, "GGAAG": 2.80, "AATAC": 2.80, "AGAGA": 2.79, "CACGT": 2.75,
        "ATAAC": 2.73, "ATAAA": 2.73, "CACGA": 2.71, "AACAG": 2.68, "CGCAT": 2.67,
        "TGCAG": 2.66, "TATTC": 2.64, "TGCAA": 2.61, "CACTA": 2.59, "CACCA": 2.55,
        "AGAGT": 2.55, "TACGT": 2.55, "TCGTA": 2.55, "AGATT": 2.55, "TACTA": 2.52,
        "GCGAA": 2.52, "TCGGA": 2.49, "TCCAT": 2.48, "TACGC": 2.46, "CCGGC": 2.45,
        "GATGG": 2.44, "GGGAG": 2.42, "TTTAG": 2.36, "CCGGT": 2.35, "AGAGG": 2.29,
        "AGATC": 2.27, "AGACT": 2.26, "TGGAA": 2.23, "TACTT": 2.23, "AGACA": 2.22,
        "TCGTC": 2.19, "TACGG": 2.18, "TACCA": 2.17, "CGCAC": 2.16, "CCGTT": 2.14,
        "CCGCG": 2.12, "GATGC": 2.11, "TGGGT": 2.07, "GATTG": 2.07, "GATTT": 2.07,
        "TCGGG": 2.04, "TTCAT": 2.03, "TTTAA": 2.03, "CTCAC": 2.00, "GATTA": 2.00,
        "CTCAG": 2.00, "CACGC": 1.97, "TCGGC": 1.97, "CTGGT": 1.96, "CCCAC": 1.94,
        "CACGG": 1.93, "CCGAA": 1.92, "TCGGT": 1.92, "CACTT": 1.92, "TACTG": 1.90,
        "TCGTT": 1.89, "TGGTC": 1.88, "AGGAC": 1.88, "CCGCA": 1.88, "TTTAC": 1.86,
        "AGACG": 1.86, "CCTAA": 1.84, "CCGCT": 1.83, "CACTG": 1.82, "AGGAT": 1.81,
        "CCCAG": 1.80, "CGGAA": 1.80, "GATGA": 1.80, "GGTAA": 1.79, "TACCG": 1.78,
        "TCCAC": 1.76, "CCCAA": 1.73, "CCGTG": 1.72, "GATTC": 1.71, "GATGT": 1.70,
        "CCGTC": 1.70, "CACCG": 1.69, "CGCAG": 1.67, "AACAC": 1.66, "CCGCC": 1.62,
        "CCGTA": 1.61, "TCCAG": 1.61, "GCTAT": 1.61, "GACTA": 1.60, "AGACC": 1.58,
        "CTGCT": 1.58, "TACCC": 1.53, "CACCC": 1.53, "AGAGC": 1.50, "CTCAA": 1.47,
        "GGCAC": 1.44, "TACCT": 1.42, "GCTAG": 1.41, "GACGT": 1.41, "TCGCA": 1.40,
        "AACAA": 1.38, "CTGTT": 1.38, "TCGTG": 1.37, "CACTC": 1.34, "TGGTA": 1.33,
        "GGGAA": 1.33, "TACTC": 1.33, "CTGTG": 1.31, "GACGC": 1.30, "TCGCT": 1.30,
        "GACGA": 1.29, "CTGGC": 1.27, "TCCAA": 1.26, "CGCAA": 1.22, "AGATG": 1.22,
        "CGGGT": 1.21, "TCGCC": 1.21, "CTGTA": 1.19, "GGCAT": 1.18, "GACGG": 1.15,
        "GTGAA": 1.13, "TGGGC": 1.10, "CGGTC": 1.08, "CGGGC": 1.07, "CTGAA": 1.05,
        "GGCAG": 1.03, "TGGTT": 0.99, "CACCT": 0.99, "CTGCA": 0.98, "GACCA": 0.97,
        "TCGCG": 0.95, "CTGCG": 0.94, "GCTAA": 0.93, "GACCG": 0.89, "CTGTC": 0.87,
        "TGGCC": 0.85, "CTGCC": 0.85, "TTCAC": 0.84, "GCTAC": 0.77, "GACCT": 0.75,
        "TTCAA": 0.73, "GACTG": 0.72, "TGGTG": 0.71, "CGGTT": 0.68, "TTCAG": 0.66,
        "CCGGG": 0.60, "GGGTC": 0.58, "CGGTG": 0.55, "TTGGT": 0.55, "GACCC": 0.54,
        "GGCAA": 0.50, "CGGTA": 0.49, "TTGTA": 0.47, "TGGCA": 0.46, "GTTAT": 0.44,
        "GGGGT": 0.41, "TTGGC": 0.37, "CTGGG": 0.37, "GACTC": 0.37, "GGGGC": 0.37,
        "CGGCC": 0.35, "TTGTT": 0.32, "TTGGG": 0.30, "GGGTG": 0.27, "AAGGT": 0.26,
        "TTGTC": 0.22, "AGGAG": 0.21, "GGGCA": 0.20, "TGGGG": 0.17, "GCCAT": 0.15,
        "CGGCA": 0.13, "CCGGA": 0.11, "TTGCT": 0.08, "GTCAT": 0.06, "TGGCG": 0.06,
        "AGTAT": 0.05, "TTGCG": 0.05, "TTGTG": 0.04, "CGGGG": 0.01, "AAGCT": 0.00,
        "TTGGA": -0.01, "AAGTA": -0.01, "TTGCC": -0.02, "GGGTA": -0.05, "TTGCA": -0.09,
        "GGGCC": -0.11, "GTTAG": -0.13, "CGGCG": -0.17, "GCGTA": -0.21, "CTGGA": -0.26,
        "AAGCG": -0.27, "AAGTC": -0.28, "GCGGT": -0.31, "AAGTG": -0.33, "GTTAC": -0.33,
        "GCCAA": -0.34, "GCGGG": -0.35, "GGGGG": -0.36, "AAGCA": -0.39, "CCTCA": -0.39,
        "GCCAC": -0.40, "GCGTT": -0.40, "GCGTG": -0.41, "GTTAA": -0.42, "TGGGA": -0.43,
        "AAGCC": -0.44, "AAGTT": -0.44, "TGTCA": -0.44, "AAGGG": -0.45, "CCTCG": -0.45,
        "CCTCT": -0.45, "TGTCT": -0.46, "CCTTT": -0.48, "ACAGT": -0.48, "GTCAC": -0.49,
        "ACAGA": -0.49, "CGGGA": -0.49, "GCCAG": -0.50, "AGTAG": -0.50, "GCGCG": -0.51,
        "AGGAA": -0.51, "GCGGA": -0.51, "GGGCG": -0.52, "GGGTT": -0.53, "ACGAT": -0.54,
        "AAGGA": -0.54, "GTGGT": -0.55, "CCTTA": -0.56, "GCGGC": -0.56, "GCGCA": -0.56,
        "GTCAG": -0.56, "AAGGC": -0.56, "CCTCC": -0.57, "TCTTA": -0.57, "CCTGT": -0.59,
        "TGTGC": -0.59, "TGGCT": -0.59, "AGTAC": -0.59, "TCTGC": -0.59, "ACGAG": -0.60,
        "GGTCG": -0.60, "GCGTC": -0.60, "ACAGG": -0.61, "CCTTC": -0.61, "ACATA": -0.62,
        "GGTCC": -0.62, "TCTTT": -0.62, "TGTCC": -0.62, "TCTTC": -0.63, "GGTCA": -0.63,
        "AGCAT": -0.64, "TCTGT": -0.64, "TGTGT": -0.64, "GCGCT": -0.65, "TCTCG": -0.65,
        "TGTCG": -0.66, "TCTCC": -0.67, "GACTT": -0.68, "TCTCA": -0.68, "CCTTG": -0.68,
        "GGGCT": -0.69, "TCTTG": -0.69, "GGGGA": -0.69, "TGTGG": -0.69, "CTTCG": -0.69,
        "ACATG": -0.69, "ACGAC": -0.70, "ATGAT": -0.71, "CGGCT": -0.72, "CGTGG": -0.72,
        "TGTGA": -0.73, "GCGCC": -0.74, "ACATT": -0.74, "TGCGA": -0.74, "GTGGC": -0.74,
        "TCTCT": -0.75, "CTTCA": -0.75, "GTGGG": -0.76, "GTGTG": -0.76, "CGCCA": -0.76,
        "ACACT": -0.77, "AGCAC": -0.77, "CGCGG": -0.77, "CTTCC": -0.77, "TGTTT": -0.77,
        "ACACC": -0.78, "CCTGC": -0.78, "TGTTC": -0.78, "CTTCT": -0.78, "AGTAA": -0.79,
        "CCCGC": -0.79, "ACACA": -0.80, "CGTCC": -0.80, "GTGTT": -0.80, "CGTCT": -0.81,
        "GGTGT": -0.82, "ACAGC": -0.82, "ACGAA": -0.82, "GCCGC": -0.83, "CGTGC": -0.83,
        "AGCAA": -0.83, "TGTTG": -0.83, "CTTTT": -0.83, "CCCGT": -0.84, "TGTTA": -0.84,
        "TCCCC": -0.84, "GTCAA": -0.84, "CCCTA": -0.84, "GTGCG": -0.84, "AATCT": -0.84,
        "CTTGT": -0.85, "GTGCT": -0.85, "CTTGC": -0.85, "GGTCT": -0.86, "AGCAG": -0.86,
        "TCCGT": -0.86, "CTCGC": -0.87, "TGCTA": -0.87, "ACATC": -0.87, "ATATA": -0.87,
        "CCCCC": -0.87, "GTGCA": -0.87, "AGGCA": -0.87, "GTGTC": -0.88, "CGTCA": -0.88,
        "CGTTC": -0.88, "TGCCA": -0.88, "CCCCG": -0.88, "TCCTT": -0.89, "TGCGC": -0.89,
        "AATGG": -0.89, "CTCGA": -0.89, "GTGCC": -0.89, "CCCGG": -0.89, "TCCGC": -0.89,
        "ACGGT": -0.89, "TCCTA": -0.89, "CGTGA": -0.89, "CTTTG": -0.89, "ATATT": -0.90,
        "TCCGG": -0.90, "ATAGA": -0.90, "CGTCG": -0.90, "ACCGC": -0.90, "ATAGG": -0.90,
        "TTCGA": -0.90, "ATACC": -0.90, "TGCCG": -0.90, "CCCTC": -0.90, "GTGGA": -0.90,
        "GGTGG": -0.91, "TTTCG": -0.91, "TGCGT": -0.91, "TTCGG": -0.91, "TGCTG": -0.91,
        "CTCTA": -0.91, "TCCCG": -0.91, "CGTGT": -0.91, "CTCTT": -0.92, "TGCCC": -0.92,
        "TCCCT": -0.92, "CGCTC": -0.92, "ATAGT": -0.92, "GGCGA": -0.92, "ATGAC": -0.92,
        "ACACG": -0.92, "AATCA": -0.92, "TCCGA": -0.93, "GGCGC": -0.93, "TCTGG": -0.93,
        "CCCTT": -0.93, "AATCG": -0.93, "AGGCT": -0.93, "GGTGC": -0.93, "TTTCT": -0.94,
        "AACTC": -0.94, "ACCTA": -0.94, "TTTCA": -0.94, "GCTGG": -0.94, "CTTTA": -0.94,
        "CTCCG": -0.94, "ATAGC": -0.94, "TTCCC": -0.94, "ATACG": -0.94, "AATGC": -0.94,
        "TTTGC": -0.94, "CGCCT": -0.95, "AGTCC": -0.95, "CGTTT": -0.95, "ATGAG": -0.95,
        "ACTAC": -0.95, "CCCGA": -0.95, "GGTTA": -0.96, "TGCCT": -0.96, "TTCCG": -0.96,
        "TCTGA": -0.96, "CGCGA": -0.96, "ATGCG": -0.96, "GCTCC": -0.97, "TCCTG": -0.97,
        "ACGCA": -0.97, "AGTGC": -0.97, "GCTTT": -0.97, "CGTTG": -0.97, "TTTGG": -0.97,
        "TTTGT": -0.97, "TGCTC": -0.97, "GCTCA": -0.97, "GGCGG": -0.98, "ACTGT": -0.98,
        "GTTGG": -0.98, "AATTA": -0.98, "AATCC": -0.98, "GGCCA": -0.98, "TGCTT": -0.98,
        "CTCGT": -0.98, "GTGTA": -0.98, "ACCAT": -0.98, "CTCCT": -0.98, "GCCCC": -0.98,
        "AGGCC": -0.99, "AATTG": -0.99, "CCCCA": -0.99, "AATGT": -0.99, "CCTGA": -0.99,
        "TTCCT": -0.99, "CTTGG": -0.99, "GCCCA": -0.99, "ACTTA": -0.99, "ACTAT": -1.00,
        "CTTTC": -1.00, "TTCCA": -1.00, "ACGCC": -1.00, "TTTGA": -1.00, "ACCAA": -1.00,
        "GGCTA": -1.00, "AGGGA": -1.00, "GCTTA": -1.00, "AATTC": -1.00, "AGCGC": -1.00,
        "ACTAG": -1.01, "TCCTC": -1.01, "GGTGA": -1.01, "TTTTA": -1.01, "AGGGG": -1.01,
        "CGTTA": -1.01, "AATTT": -1.01, "TTCTG": -1.01, "AGTGG": -1.01, "ATTTG": -1.01,
        "AATGA": -1.01, "TTCTA": -1.01, "GCCCG": -1.01, "ATACA": -1.01, "ACCCG": -1.01,
        "ACCTC": -1.01, "GTCCC": -1.01, "TGCGG": -1.01, "ACCAG": -1.01, "GGTTG": -1.02,
        "AGTCA": -1.02, "CTCTC": -1.02, "ATGGG": -1.02, "ATGCT": -1.02, "CCCTG": -1.02,
        "GCTCT": -1.02, "AGGCG": -1.02, "ATGGC": -1.02, "CTCTG": -1.02, "ACGTA": -1.02,
        "ATGGT": -1.02, "ACTCT": -1.03, "ATTCG": -1.03, "GCCTT": -1.03, "GTCCG": -1.03,
        "TTCTT": -1.03, "ACCTG": -1.03, "GGTTT": -1.03, "ACTGA": -1.03, "AGGTG": -1.03,
        "ATGAA": -1.03, "ATCGG": -1.03, "ATTCC": -1.03, "ACCCA": -1.03, "ACTCA": -1.03,
        "GTTTG": -1.03, "CGCTA": -1.03, "CTCCC": -1.03, "AACTG": -1.03, "ATCAG": -1.03,
        "GCTTG": -1.04, "TCCCA": -1.04, "AACCG": -1.04, "AACCC": -1.04, "GCTGT": -1.04,
        "TTCGT": -1.04, "ACCCT": -1.04, "AGCCG": -1.04, "TTTCC": -1.04, "ATATG": -1.04,
        "ACTTT": -1.04, "AACGT": -1.04, "ATGGA": -1.04, "GTCCT": -1.04, "GTTGA": -1.04,
        "ATTGT": -1.05, "AGTGA": -1.05, "ATCGT": -1.05, "ATCGC": -1.05, "ACTCG": -1.05,
        "GCCTG": -1.05, "ATTAC": -1.05, "ACCAC": -1.05, "GCTTC": -1.05, "TTTTT": -1.05,
        "AACGG": -1.05, "CTTGA": -1.05, "ATCAT": -1.05, "GGCCT": -1.05, "GTCGG": -1.05,
        "GTCCA": -1.05, "ATTAT": -1.06, "CTCCA": -1.06, "GTCTG": -1.06, "AACCT": -1.06,
        "ATTGG": -1.06, "TTTTC": -1.06, "TTCGC": -1.06, "ATGCC": -1.06, "GTTCT": -1.06,
        "ACTGC": -1.06, "GGCCG": -1.06, "ATATC": -1.06, "CCTGG": -1.06, "AGTTA": -1.06,
        "GCTGA": -1.06, "ATCCT": -1.06, "GCTGC": -1.06, "GGCTT": -1.07, "AACCA": -1.07,
        "GTTGT": -1.07, "AGCCC": -1.07, "AACTT": -1.07, "ACGTT": -1.07, "AGGTA": -1.07,
        "ACGTC": -1.07, "ACTCC": -1.07, "GTCGC": -1.07, "AGTCG": -1.07, "GTCGA": -1.07,
        "AACGA": -1.07, "ATCTT": -1.07, "ATTAG": -1.07, "AGGGT": -1.08, "ATGTA": -1.08,
        "AGTTT": -1.08, "ATCTC": -1.08, "AGCCA": -1.08, "ATGTG": -1.08, "CGCCC": -1.08,
        "AGGGC": -1.08, "GTCGT": -1.08, "GTTTT": -1.08, "AACTA": -1.08, "GGCGT": -1.08,
        "GGTTC": -1.09, "ACGGG": -1.09, "GTTCG": -1.09, "GTTTC": -1.09, "ATCCC": -1.09,
        "ATTAA": -1.09, "GCTCG": -1.09, "ACTTC": -1.09, "TTTTG": -1.09, "ATTCT": -1.09,
        "ATTGA": -1.10, "GTTTA": -1.10, "ATCTA": -1.10, "AGCTG": -1.10, "ATCGA": -1.10,
        "AGCTC": -1.10, "GCCGT": -1.10, "ACGGC": -1.10, "ACTTG": -1.10, "CGCGT": -1.10,
        "GTTCC": -1.10, "ATGTT": -1.11, "ATACT": -1.11, "GCCGG": -1.11, "ACGCG": -1.11,
        "GGCCC": -1.11, "ACGTG": -1.11, "AGTTG": -1.11, "ATGTC": -1.11, "AGGTT": -1.11,
        "GGCTG": -1.11, "CCCCT": -1.11, "GCCTC": -1.11, "GTTGC": -1.12, "ATTCA": -1.12,
        "AGTGT": -1.12, "ATTTC": -1.12, "ATTTT": -1.12, "GTCTC": -1.12, "GCCGA": -1.12,
        "ATCCA": -1.12, "CTCGG": -1.12, "ACCGT": -1.13, "CGCTT": -1.13, "ACTAA": -1.13,
        "ATCCG": -1.13, "ATTGC": -1.13, "ACGCT": -1.13, "AGGTC": -1.13, "ACTGG": -1.13,
        "AACGC": -1.13, "GCCCT": -1.14, "ACCGG": -1.14, "ATCAC": -1.14, "GTTCA": -1.14,
        "AGCTT": -1.14, "AGTTC": -1.15, "ACGGA": -1.15, "GCCTA": -1.15, "ATCAA": -1.15,
        "ATCTG": -1.15, "ACCTT": -1.16, "ACCCC": -1.16, "AGCGT": -1.16, "TTCTC": -1.16,
        "AGTCT": -1.17, "ATGCA": -1.18, "GTCTA": -1.19, "CGCTG": -1.19, "ATTTA": -1.19,
        "ACCGA": -1.19, "AGCGA": -1.21, "GTCTT": -1.21, "AGCTA": -1.21, "CGCCG": -1.22,
        "AGCGG": -1.22, "GGCTC": -1.23, "AGCCT": -1.23, "CGCGC": -1.25
    }
    return scoring_table.get(kmer, 0.0)


def get_alternative_pfs_score(kmer):
    """Returns the alternative PFS score using the positional matrix."""
    scoring_table = {
        0: {'A': 0.6, 'C': 1.2, 'G': 0.99, 'T': 1.2},
        1: {'A': 2, 'C': 0.86, 'G': 0.87, 'T': 0.84},
        2: {'A': 3.7, 'C': 0.7, 'G': 1.0, 'T': 0.8},
        3: {'A': 1.9, 'C': 0.83, 'G': 0.85, 'T': 0.85},
        4: {'A': 1, 'C': 1, 'G': 1, 'T': 1},
    }
    score = 1.0
    for i, base in enumerate(kmer):
        base_score = scoring_table[i].get(base)
        if base_score is None:
            return 0.0  # Return 0 gracefully if bad characters exist
        score *= base_score
    return score


def score_kmer(kmer, method="published"):
    kmer = kmer.upper()

    if len(kmer) < 4:
        return 0.0
    elif len(kmer) > 5:
        return 0.0  # Fallback gracefully instead of ValueError to support batch runs

    pub_score = get_published_pfs_score(kmer)
    alt_score = get_alternative_pfs_score(kmer)

    if method == "published":
        return round(pub_score, 1)
    elif method == "alternative":
        return round(alt_score, 1)
    elif method == "average":
        return round((pub_score + alt_score) / 2, 1)
    else:
        return round(pub_score, 1)


def analyse_linear_sequence(sequence, terms=None, gblock_tail=None, pfs_method="published"):
    sequence = sequence.upper()
    terms = terms if terms is not None else DEFAULT_TERMS
    gblock_tail = gblock_tail if gblock_tail is not None else DEFAULT_GBLOCK_TAIL

    results = []

    for term in terms:
        for start_index in [i for i in range(len(sequence)) if sequence.startswith(term, i)]:
            upstream_start = max(0, start_index - 28)
            target_region = sequence[upstream_start:start_index]

            if len(target_region) != 28:
                continue

            # Reverse complement
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
            rc_target = "".join(complement.get(base, base) for base in reversed(target_region))

            gc_percent = round(GC_calc(rc_target), 1)

            adjusted_term = term
            if len(adjusted_term) < 5:
                pad_right = sequence[start_index + len(adjusted_term):start_index + 5]
                adjusted_term += pad_right

            pfs_score = score_kmer(adjusted_term, method=pfs_method)

            # Construct gBlock sequence with custom tail if provided
            gblock_sequence = target_region + gblock_tail

            results.append({
                "Term": term,
                "Adjusted Term": adjusted_term,
                "PFS Score": pfs_score,
                "Target Region": target_region,
                "RC Target Region": rc_target,
                "GC % in sequence": gc_percent,
                "gBlock": gblock_sequence
            })

    df = pd.DataFrame(results)

    if df.empty:
        return df

    # Sort top-down by PFS score and drop duplicates based on the Target Region
    df = df.sort_values(by="PFS Score", ascending=False)
    df = df.drop_duplicates(subset=["Target Region"], keep="first")

    return df


def analyse_sequence(upstream, downstream, terms=None, linearoverhang="",
                     up_length=27, down_length=30, gblock_tail=None,
                     pfs_method="published", min_dist=-7, max_dist=17):
    terms = terms if terms is not None else DEFAULT_TERMS
    gblock_tail = gblock_tail if gblock_tail is not None else DEFAULT_GBLOCK_TAIL

    length_downstream = len(downstream)
    joint_sequence = downstream + upstream
    truncated_sequence = joint_sequence[length_downstream - down_length: length_downstream + up_length]

    truncated_upstream = upstream[:up_length]
    truncated_downstream = downstream[-17:]
    truncated_linear_overhang = linearoverhang[-28:]

    bsj_site = downstream[-1] + upstream[0]
    search_sequence = downstream[-4:] + upstream[:up_length]
    length_diff = len(joint_sequence) - len(upstream)

    results = []

    for term in terms:
        for start_index_trunc in [i for i in range(len(search_sequence)) if search_sequence.startswith(term, i)]:
            bsj_distance = start_index_trunc - 4

            # Apply distance filtering constraints (exclusive to Circular RNA mode)
            if not (min_dist <= bsj_distance <= max_dist):
                continue

            within_bsj = bsj_distance <= 1
            start_index_long = start_index_trunc + length_diff

            target_start = max(0, start_index_long - 28)
            target_sequence = joint_sequence[target_start:start_index_long - 4]

            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
            rc_target_sequence = "".join(complement.get(base, base) for base in reversed(target_sequence))
            gc_percent = round(GC_calc(rc_target_sequence), 1)

            overhang_sequence = joint_sequence[target_start:length_downstream]
            linear_overlap_percent = round(compare_sequence_similarity(truncated_linear_overhang, overhang_sequence), 1)

            adjusted_sequence = term
            if len(adjusted_sequence) < 5:
                pad_right = search_sequence[start_index_trunc + len(adjusted_sequence):start_index_trunc + 5]
                adjusted_sequence += pad_right

            pfs_score = score_kmer(adjusted_sequence, method=pfs_method)

            # Use custom gBlock tail if provided
            gblock_sequence = target_sequence + gblock_tail

            results.append({
                "Term": term,
                "Adjusted Term": adjusted_sequence,
                "PFS Score": pfs_score,
                "Within BSJ": within_bsj,
                "Closest distance from BSJ": bsj_distance,
                "Target Region": target_sequence,
                "RC Target Region": rc_target_sequence,
                "GC % in sequence": gc_percent,
                "Linear Sequence Overhang Similarity Percent": linear_overlap_percent,
                "gBlock": gblock_sequence
            })

    df = pd.DataFrame(results)

    if df.empty:
        return truncated_sequence, truncated_upstream, truncated_downstream, bsj_site, df

    # Sort top-down by PFS score, then drop duplicate distances (keeping the highest score)
    df = df.sort_values(by="PFS Score", ascending=False)
    df = df.drop_duplicates(subset=["Closest distance from BSJ"], keep="first")

    return truncated_sequence, truncated_upstream, truncated_downstream, bsj_site, df


def prompt_pfs_method():
    print("\nSelect PFS Scoring Method:")
    print("1. Published Data")
    print("2. Alternative Method")
    print("3. Average of Both")
    choice = input("Enter choice (1-3): ").strip()

    if choice == '2':
        return "alternative"
    elif choice == '3':
        return "average"
    return "published"


def prompt_circ_bounds():
    print("\nSelect BSJ Distance Filtering Bounds:")
    print("1. Default (-7 to +17)")
    print("2. Custom Limits")
    choice = input("Enter choice (1-2): ").strip()

    if choice == '2':
        try:
            min_dist = int(input("Enter minimum distance (e.g., -7): ").strip())
            max_dist = int(input("Enter maximum distance (e.g., 17): ").strip())
            return min_dist, max_dist
        except ValueError:
            print("Invalid input. Using defaults.")

    return -7, 17


def main():
    terms = DEFAULT_TERMS
    gblock_tail = DEFAULT_GBLOCK_TAIL

    while True:
        print("\n=== Sequence Analyzer ===")
        print("1. Circular RNA Mode")
        print("2. Linear Sequence Mode")
        print("3. Custom Terms Mode")
        mode = input("Select analysis mode (1, 2, or 3): ").strip()

        if mode == '1':
            print("\n=== Circular RNA Sequence Analyzer ===")
            pfs_method = prompt_pfs_method()
            min_dist, max_dist = prompt_circ_bounds()

            upstream = input("\nEnter upstream sequence: ").strip().upper()
            downstream = input("Enter downstream sequence: ").strip().upper()
            linearoverhang = input("Enter upstream of upstream sequence: ").strip().upper()

            truncated_seq, truncated_upstream, truncated_downstream, bsj_site, df_results = analyse_sequence(
                upstream, downstream, terms, linearoverhang,
                pfs_method=pfs_method, min_dist=min_dist, max_dist=max_dist)

            print(f"\nTruncated Sequence: {truncated_seq}")
            print(f"Truncated Downstream: {truncated_downstream}")
            print(f"Truncated Upstream: {truncated_upstream}")
            print(f"BSJ Site: {bsj_site}\n")

        elif mode == '2':
            print("\n=== Linear Sequence Analyzer ===")
            pfs_method = prompt_pfs_method()

            linear_seq = input("\nEnter the linear DNA sequence: ").strip().upper()
            df_results = analyse_linear_sequence(linear_seq, terms, pfs_method=pfs_method)

        elif mode == '3':
            custom_input = input("\nEnter your custom terms separated by commas (e.g. TAGAC,TGACA): ").strip().upper()
            custom_terms = [term.strip() for term in custom_input.split(",") if term.strip()]

            if not custom_terms:
                print("No valid terms entered. Returning to main menu.\n")
                continue

            print("\nChoose analysis mode for custom terms:")
            print("1. Circular RNA Mode")
            print("2. Linear Sequence Mode")
            custom_mode = input("Select (1 or 2): ").strip()

            if custom_mode == '1':
                print("\n=== Circular RNA Sequence Analyzer (Custom Terms) ===")
                pfs_method = prompt_pfs_method()
                min_dist, max_dist = prompt_circ_bounds()

                upstream = input("\nEnter upstream sequence: ").strip().upper()
                downstream = input("Enter downstream sequence: ").strip().upper()
                linearoverhang = input("Enter upstream of upstream sequence: ").strip().upper()

                truncated_seq, truncated_upstream, truncated_downstream, bsj_site, df_results = analyse_sequence(
                    upstream, downstream, custom_terms, linearoverhang,
                    pfs_method=pfs_method, min_dist=min_dist, max_dist=max_dist)

                print(f"\nTruncated Sequence: {truncated_seq}")
                print(f"Truncated Downstream: {truncated_downstream}")
                print(f"Truncated Upstream: {truncated_upstream}")
                print(f"BSJ Site: {bsj_site}\n")

            elif custom_mode == '2':
                print("\n=== Linear Sequence Analyzer (Custom Terms) ===")
                pfs_method = prompt_pfs_method()

                linear_seq = input("\nEnter the linear DNA sequence: ").strip().upper()
                df_results = analyse_linear_sequence(linear_seq, custom_terms, pfs_method=pfs_method)

            else:
                print("Invalid option for custom analysis mode. Returning to main menu.\n")
                continue

        else:
            print("Invalid option. Please enter 1, 2, or 3.\n")
            continue

        print("Search Results:\n")
        if df_results.empty:
            print("No matching motifs found in the sequence.")
        else:
            print(df_results.to_string(index=False))

            # CSV Output hook
            export = input("\nWould you like to save these results to a CSV? (y/n): ").strip().lower()
            if export == 'y':
                filename = input("Enter a filename (e.g., results.csv): ").strip()
                if not filename.endswith('.csv'):
                    filename += '.csv'
                df_results.to_csv(filename, index=False)
                print(f"Results saved successfully to {filename}.")

        restart = input("\nWould you like to run another sequence? (y/n): ").strip().lower()
        if restart != 'y':
            print("\nExiting the program. Goodbye!")
            break


if __name__ == "__main__":
    main()