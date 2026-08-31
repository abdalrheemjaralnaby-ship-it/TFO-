"""
specificity.py
Checks exact binding uniqueness inside the selected target region.
"""
from typing import List, Dict


def find_exact_occurrences(target_sequence: str, full_sequence: str) -> Dict:
    """
    Count exact occurrences of a candidate target sequence in the selected region.
    Uniqueness is true only when the candidate binds at one exact position.
    """
    target_upper = target_sequence.upper()
    seq_upper = full_sequence.upper()
    target_len = len(target_upper)
    positions: List[int] = []

    if not target_upper or target_len > len(seq_upper):
        return {"match_count": 0, "positions": [], "is_unique": False}

    for i in range(len(seq_upper) - target_len + 1):
        if seq_upper[i:i + target_len] == target_upper:
            positions.append(i)

    return {
        "match_count": len(positions),
        "positions": positions,
        "is_unique": len(positions) == 1,
    }
