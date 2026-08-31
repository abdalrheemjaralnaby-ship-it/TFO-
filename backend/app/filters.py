"""
filters.py
Applies dynamic filters to a list of TFOCandidate dicts.
"""
from typing import List, Dict
from .models import FilterOptions


def apply_filters(candidates: List[Dict], filters: FilterOptions) -> List[Dict]:
    """
    Filter generated TFO candidates with the core biological filters:
    target mismatches, TFO purine content, then exact unique binding.
    """
    result = []
    
    for c in candidates:
        if c["target_mismatches"] > filters.max_target_mismatches:
            continue

        if c["purine_count"] > filters.max_tfo_purines:
            continue

        if filters.only_unique and not c["is_unique"]:
            continue

        result.append(c)

    return result
