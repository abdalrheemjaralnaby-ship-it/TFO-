"""
tfo_generator.py
Generates multiple TFO candidates per TTS using sliding windows of
short / medium / long lengths.
"""
from typing import List, Optional
from .models import TTSResult

# Length ranges per category
LENGTH_RANGES = {
    "short":  range(12, 16),   # 12–15 nt
    "medium": range(16, 26),   # 16–25 nt
    "long":   range(26, 41),   # 26–40 nt
}

COMPLEMENT = str.maketrans("ATGCN", "TACGN")


def _complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)


def _make_tfo(tts_segment: str) -> str:
    """
    Convert a TTS segment to its TFO sequence.
    """
    return _complement(tts_segment)


def generate_tfos(
    tts_list: List[TTSResult],
    length_categories: Optional[List[str]] = None,
) -> List[dict]:
    """
    For each TTS region, produce overlapping window TFO candidates.
    Returns a list of raw candidate dicts.
    """
    if length_categories is None:
        length_categories = ["short", "medium", "long"]

    candidates = []

    for tts in tts_list:
        tts_seq = tts.sequence
        tts_len = len(tts_seq)
        local_count = 0

        for category in length_categories:
            for win_len in LENGTH_RANGES.get(category, []):
                if win_len > tts_len:
                    continue
                # Slide across the TTS with step=1 (fully overlapping)
                for offset in range(0, tts_len - win_len + 1):
                    segment = tts_seq[offset: offset + win_len]
                    tfo_seq = _make_tfo(segment)
                    abs_start = tts.start + offset
                    abs_end = abs_start + win_len - 1
                    
                    candidates.append({
                        "tts_start": tts.start,
                        "tts_end": tts.end,
                        "tts_sequence": segment,  # Match the candidate window
                        "tfo_sequence": tfo_seq,
                        "start": abs_start,
                        "end": abs_end,
                        "length": win_len,
                        "length_category": category,
                    })
                    local_count += 1
        
        # Robustness: ensure each TTS produces at least one candidate if it made it here
        if local_count == 0 and tts_len > 0:
            tfo_seq = _make_tfo(tts_seq)
            candidates.append({
                "tts_start": tts.start,
                "tts_end": tts.end,
                "tts_sequence": tts_seq,
                "tfo_sequence": tfo_seq,
                "start": tts.start,
                "end": tts.end,
                "length": tts_len,
                "length_category": "short" if tts_len < 16 else "medium",
            })

    return candidates
