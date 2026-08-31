from pydantic import BaseModel, Field
from typing import Optional, List


# ─── Input Models ────────────────────────────────────────────────────────────

class ParseRequest(BaseModel):
    input_type: str  # "manual" | "fasta"
    sequence: Optional[str] = None
    fasta_content: Optional[str] = None


class RegionSelection(BaseModel):
    start: int
    end: int
    label: Optional[str] = "custom"
    feature_type: Optional[str] = None


class FilterOptions(BaseModel):
    only_unique: bool = False
    max_target_mismatches: int = 2
    max_tfo_purines: int = 100
    length_categories: List[str] = ["short", "medium", "long"]


class TFOFindRequest(BaseModel):
    sequence: Optional[str] = ""
    region: Optional[RegionSelection] = None
    region_type: Optional[str] = None  # "full", "manual", "annotation", "promoter"
    promoter_length: int = 1000
    tts_min_length: int = 15
    tts_purine_ratio: float = 0.8
    filters: Optional[FilterOptions] = None
    is_ncbi: bool = False
    accession: Optional[str] = None


# ─── Output Models ───────────────────────────────────────────────────────────

class AnnotationFeature(BaseModel):
    type: str
    label: str
    start: int
    end: int
    strand: Optional[int] = 1


class ParseResponse(BaseModel):
    sequence: str
    length: int
    annotations: List[AnnotationFeature] = []
    source: str = "manual"


class TTSResult(BaseModel):
    start: int
    end: int
    sequence: str
    length: int
    purine_ratio: float


class TFOCandidate(BaseModel):
    tts_start: int
    tts_end: int
    tts_sequence: str
    tfo_sequence: str
    start: int
    end: int
    length: int
    match_count: int
    positions: List[int]
    is_unique: bool
    target_mismatches: int
    purine_count: int = 0
    gc_content: float
    score: float
    length_category: str


class TFOFindResponse(BaseModel):
    tts_regions: List[TTSResult]
    candidates: List[TFOCandidate]
    total_tts: int
    total_candidates: int
    unique_candidates: int
    sequence_length: int
    target_sequence: str = ""
    target_offset: int = 0
