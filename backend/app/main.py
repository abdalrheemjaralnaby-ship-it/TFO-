"""
main.py
FastAPI application with three endpoints:
  POST  /sequence/parse         – parse manual or FASTA input
  GET   /sequence/fetch/{id}    – fetch from NCBI by accession
  POST  /tfo/find               – full TTS → TFO pipeline
"""
from fastapi import FastAPI, HTTPException, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import RedirectResponse
from typing import Optional
import traceback

from .models import (
    ParseRequest, ParseResponse,
    TFOFindRequest, TFOFindResponse,
    TFOCandidate, TTSResult,
    FilterOptions,
)
from .sequence_parser import parse_manual, parse_fasta_content, fetch_ncbi, search_ncbi
from .tts_detector import detect_tts
from .tfo_generator import generate_tfos
from .specificity import find_exact_occurrences
from .scoring import (
    count_mismatches, compute_gc, compute_tfo_score, count_tfo_purines
)
from .filters import apply_filters

app = FastAPI(
    title="TriplexLab – TFO Finder API",
    description="Backend API for detecting TTS regions and generating TFO candidates.",
    version="1.0.0",
)

# ─── CORS ─────────────────────────────────────────────────────────────────────
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ─── FRONTEND STATIC FILES ───────────────────────────────────────────────────
import os
try:
    frontend_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../frontend"))
    if os.path.exists(frontend_dir):
        # We need a dedicated root route so it catches explicitly
        @app.get("/")
        async def root():
            return RedirectResponse(url="/ui/index.html")
            
        app.mount("/ui", StaticFiles(directory=frontend_dir), name="ui")
except Exception as e:
    pass


# ─── POST /sequence/parse ─────────────────────────────────────────────────────
@app.post("/sequence/parse", response_model=ParseResponse)
async def parse_sequence(request: ParseRequest):
    """
    Parse a DNA sequence from manual text or FASTA content.
    """
    try:
        if request.input_type == "manual":
            if not request.sequence:
                raise HTTPException(status_code=400, detail="No sequence provided.")
            seq, annotations = parse_manual(request.sequence)
        elif request.input_type == "fasta":
            if not request.fasta_content:
                raise HTTPException(status_code=400, detail="No FASTA content provided.")
            seq, annotations = parse_fasta_content(request.fasta_content)
        else:
            raise HTTPException(status_code=400, detail=f"Unknown input_type: {request.input_type}")

        return ParseResponse(
            sequence=seq,
            length=len(seq),
            annotations=annotations,
            source=request.input_type,
        )
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Parse error: {traceback.format_exc()}")


# ─── GET /sequence/search ─────────────────────────────────────────────────────
@app.get("/sequence/search")
async def search_sequence(query: str, organism: str = ""):
    """
    Detects input type (NM_, NC_, or Gene name) and returns relevant transcripts.
    """
    try:
        return search_ncbi(query, organism)
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"NCBI search error: {str(e)}")


# ─── GET /sequence/fetch/{id} ─────────────────────────────────────────────────
@app.get("/sequence/fetch/{accession}", response_model=ParseResponse)
async def fetch_sequence(accession: str, types: Optional[str] = None):
    """
    Fetch a DNA sequence and annotations from NCBI by accession number or gene name.
    """
    try:
        filter_types = types.split(",") if types else None
        seq, annotations = fetch_ncbi(accession, filter_types)
        return ParseResponse(
            sequence=seq,
            length=len(seq),
            annotations=annotations,
            source="ncbi",
        )
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"NCBI fetch error: {str(e)}")


# ─── POST /tfo/find ───────────────────────────────────────────────────────────
@app.post("/tfo/find", response_model=TFOFindResponse)
async def find_tfos(request: TFOFindRequest):
    """
    Full pipeline:
      1. Extract the relevant subsequence (region or full seq)
      2. Detect TTS regions
      3. Generate TFO candidates with sliding windows
      4. Compute GC, score, and specificity for each candidate
      5. Apply filters
      6. Return results
    """
    try:
        filters = request.filters or FilterOptions()

        # --- 1. Region extraction and Validation ---
        if request.is_ncbi:
            is_promoter = (request.region_type == "promoter")

            # Promoter mode does NOT need a pre-set region; backend computes it
            if not is_promoter and not request.region:
                raise HTTPException(status_code=400, detail="Please select a region before analysis")

            if not is_promoter:
                start = request.region.start
                end = request.region.end
                if start >= end:
                    raise HTTPException(status_code=400, detail="Start must be less than end")

            if not request.accession:
                raise HTTPException(status_code=400, detail="Missing accession ID for NCBI analysis")

            try:
                from Bio import Entrez, SeqIO
                from .sequence_parser import clean_sequence, compute_promoter_region
                Entrez.email = "triplexlab@example.com"

                fetch_handle = Entrez.efetch(
                    db="nucleotide",
                    id=request.accession,
                    rettype="gb",
                    retmode="text"
                )
                record = SeqIO.read(fetch_handle, "genbank")
                fetch_handle.close()

                if is_promoter:
                    # Strategy: Try computing from the current record first.
                    # If it's a genomic record (NC_), it should work. If it's a
                    # transcript, the record commonly starts at TSS and yields 0-0.
                    promoter_seq, start, end, is_reverse, tss = compute_promoter_region(
                        record, request.promoter_length
                    )
                    working_seq = promoter_seq

                    # Fallback: resolve genomic coordinates and fetch upstream
                    # sequence from the chromosome. This is useful for NM_/XM_
                    # transcripts and other NCBI records that have no local upstream.
                    genomic_fallback_error = None
                    if not working_seq or len(working_seq) < request.promoter_length:
                        try:
                            from .sequence_parser import get_genomic_promoter
                            working_seq, start, end, is_reverse, tss = get_genomic_promoter(
                                request.accession, request.promoter_length
                            )
                        except Exception as e:
                            genomic_fallback_error = e

                    # Final fallback: analyze the available TSS-adjacent sequence
                    # from the fetched record instead of failing with a 0-0
                    # promoter. Some NCBI accessions (for example PZ_/XM_/NM_)
                    # do not expose an Entrez gene link even though the record
                    # itself has usable sequence.
                    if not working_seq:
                        seq = clean_sequence(str(record.seq))
                        fallback_len = min(len(seq), max(1, request.promoter_length))
                        if is_reverse:
                            start = max(0, len(seq) - fallback_len)
                            end = len(seq)
                            working_seq = seq[start:end]
                        else:
                            start = 0
                            end = fallback_len
                            working_seq = seq[start:end]

                    if not working_seq:
                        detail = (
                            f"Promoter analysis could not find any sequence for {request.accession}. "
                            "Try selecting a gene/CDS annotation or using a genomic accession (NC_)."
                        )
                        if genomic_fallback_error:
                            detail += f" Genomic promoter lookup failed: {genomic_fallback_error}"
                        raise HTTPException(status_code=400, detail=detail)
                else:
                    region_seq = str(record.seq)[start:end]
                    working_seq = clean_sequence(region_seq)

            except HTTPException:
                raise
            except Exception as e:
                raise HTTPException(status_code=500, detail=f"Failed to fetch sequence: {str(e)}")

            offset = start
            seq_length = len(record.seq)  # Length before slice
        else:
            if not request.sequence:
                raise HTTPException(status_code=400, detail="Sequence string missing")
            
            seq = request.sequence.upper()
            
            if request.region:
                start = max(0, request.region.start)
                end = min(len(seq), request.region.end)
                working_seq = seq[start:end]
                offset = start
            else:
                working_seq = seq
                offset = 0
            
            seq_length = len(seq)

        if not working_seq:
            raise HTTPException(status_code=400, detail="The selected region is empty.")

        # --- 2. TTS detection ---
        tts_results = detect_tts(
            working_seq,
            min_length=request.tts_min_length,
            purine_ratio_threshold=request.tts_purine_ratio,
        )

        # Adjust coordinates to absolute positions
        adjusted_tts: list[TTSResult] = []
        for tts in tts_results:
            adjusted_tts.append(TTSResult(
                start=tts.start + offset,
                end=tts.end + offset,
                sequence=tts.sequence,
                length=tts.length,
                purine_ratio=tts.purine_ratio,
            ))

        # --- 3. Generate TFO candidates ---
        raw_candidates = generate_tfos(
            adjusted_tts,
            length_categories=filters.length_categories,
        )

        # --- 4. Core candidate metrics and scoring ---
        enriched = []
        for c in raw_candidates:
            binding = find_exact_occurrences(c["tts_sequence"], working_seq)
            adjusted_positions = [pos + offset for pos in binding["positions"]]

            target_mismatches = count_mismatches(c["tfo_sequence"], c["tts_sequence"])
            if target_mismatches is None:
                continue

            tfo_purines = count_tfo_purines(c["tfo_sequence"])
            gc_content = compute_gc(c["tfo_sequence"])
            score = compute_tfo_score(c["length"], gc_content, target_mismatches)

            candidate_dict = {
                **c,
                "match_count": binding["match_count"] or 0,
                "positions": adjusted_positions,
                "is_unique": binding["is_unique"] or False,
                "target_mismatches": int(target_mismatches),
                "purine_count": int(tfo_purines),
                "gc_content": float(gc_content),
                "score": float(score)
            }
            enriched.append(candidate_dict)

        # --- 5. Apply filters ---
        filtered = apply_filters(enriched, filters)

        # --- 6. Sort by score descending (higher score = better match) ---
        filtered.sort(key=lambda x: x["score"], reverse=True)

        # Build response objects
        candidates = [TFOCandidate(**c) for c in filtered]

        return TFOFindResponse(
            tts_regions=adjusted_tts,
            candidates=candidates,
            total_tts=len(adjusted_tts),
            total_candidates=len(candidates),
            unique_candidates=sum(1 for c in candidates if c.is_unique),
            sequence_length=seq_length,
            target_sequence=working_seq,
            target_offset=offset,
        )

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"TFO pipeline error: {traceback.format_exc()}")


# ─── Health check ─────────────────────────────────────────────────────────────
@app.get("/health")
async def health():
    return {"status": "ok", "service": "TriplexLab TFO Finder API"}
