"""
tts_detector.py
Detects polypurine (TTS) regions in a DNA sequence.
Purines: A, G
"""
from typing import List
from .models import TTSResult


def _purine_ratio(segment: str) -> float:
    if not segment:
        return 0.0
    purines = sum(1 for b in segment if b in ('A', 'G'))
    return purines / len(segment)


def detect_tts(
    sequence: str,
    min_length: int = 15,
    purine_ratio_threshold: float = 0.8,
) -> List[TTSResult]:
    """
    Scan sequence with a sliding window and collect merged polypurine regions.

    Strategy:
      1. Use a window of size `min_length`.
      2. If the window meets the purine ratio, extend the region as far as
         possible while keeping the ratio above threshold.
      3. Merge overlapping / adjacent candidate regions.
    """
    seq = sequence.upper()
    n = len(seq)
    detected: List[TTSResult] = []

    if n < min_length:
        return detected

    i = 0
    while i <= n - min_length:
        window = seq[i: i + min_length]
        if _purine_ratio(window) >= purine_ratio_threshold:
            # Extend right as long as ratio stays valid
            end = i + min_length
            while end < n and _purine_ratio(seq[i: end + 1]) >= purine_ratio_threshold:
                end += 1
            region_seq = seq[i:end]
            ratio = _purine_ratio(region_seq)
            detected.append(TTSResult(
                start=i,
                end=end - 1,
                sequence=region_seq,
                length=len(region_seq),
                purine_ratio=round(ratio, 3),
            ))
            i = end  # jump past this region
        else:
            i += 1

    return _merge_overlapping(detected, seq, purine_ratio_threshold)


def _merge_overlapping(
    regions: List[TTSResult],
    seq: str,
    purine_ratio_threshold: float,
) -> List[TTSResult]:
    """Merge TTS regions that overlap or are adjacent."""
    if not regions:
        return regions
    merged: List[TTSResult] = [regions[0]]
    for current in regions[1:]:
        last = merged[-1]
        if current.start <= last.end + 1:
            new_end = max(last.end, current.end)
            new_seq = seq[last.start: new_end + 1]
            ratio = _purine_ratio(new_seq)
            merged[-1] = TTSResult(
                start=last.start,
                end=new_end,
                sequence=new_seq,
                length=len(new_seq),
                purine_ratio=round(ratio, 3),
            )
        else:
            merged.append(current)
    return merged
