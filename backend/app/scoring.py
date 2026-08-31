"""
scoring.py
Computes target mismatches, TFO purine count, GC content, and quality score.
"""


def count_mismatches(tfo: str, target: str):
    """
    Count TFO-to-target mismatches using the supported pairing rules:
    T pairs with A, and C pairs with G.
    """
    tfo = tfo.upper()
    target = target.upper()
    if len(tfo) != len(target):
        return None

    mismatches = 0
    for t, d in zip(tfo, target):
        if (t == "T" and d == "A") or (t == "C" and d == "G"):
            continue
        mismatches += 1
    return mismatches


def count_tfo_purines(tfo_sequence: str) -> int:
    """
    Count number of purine bases (A and G) directly from the TFO sequence.
    Do not transform or reverse the sequence before counting.
    """
    return tfo_sequence.upper().count("A") + tfo_sequence.upper().count("G")


def compute_gc(seq: str) -> float:
    """Return GC content as a fraction (0.0 - 1.0). Formula: (G + C) / length"""
    if not seq:
        return 0.0
    gc = sum(1 for b in seq.upper() if b in ('G', 'C'))
    return round(gc / len(seq), 4)


def compute_tfo_score(length: int, gc_content: float, target_mismatches: int) -> float:
    """
    TFO quality score based on user requirements:
    score = (length * 0.5) + (gc_content * 10) - (target_mismatches * 2)
    """
    score = (float(length) * 0.5) + (float(gc_content) * 10.0) - (float(target_mismatches) * 2.0)
    return round(max(0.0, score), 2)
