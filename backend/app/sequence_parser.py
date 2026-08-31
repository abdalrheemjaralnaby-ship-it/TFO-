"""
sequence_parser.py
Handles manual input, FASTA parsing, and NCBI GenBank fetching.
"""
import re
import io
from typing import List, Optional, Tuple
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord

from .models import AnnotationFeature

# Set a default email for NCBI Entrez (required by NCBI policy)
Entrez.email = "triplexlab@example.com"

VALID_BASES = re.compile(r'^[ACTGNacgt\s]+$')
CLEAN_RE = re.compile(r'[^ACTGNacgt]')


def clean_sequence(seq: str) -> str:
    """Uppercase and strip non-nucleotide characters."""
    return CLEAN_RE.sub('', seq.upper())


def validate_sequence(seq: str) -> Tuple[bool, str]:
    """Returns (is_valid, error_message)."""
    if not seq:
        return False, "Sequence is empty."
    cleaned = clean_sequence(seq)
    if len(cleaned) == 0:
        return False, "No valid nucleotide characters found."
    invalid = set(cleaned) - {'A', 'C', 'T', 'G', 'N'}
    if invalid:
        return False, f"Invalid characters after cleaning: {invalid}"
    return True, ""


def parse_manual(raw: str) -> Tuple[str, List[AnnotationFeature]]:
    """Parse and validate a manually entered sequence."""
    cleaned = clean_sequence(raw)
    is_valid, err = validate_sequence(cleaned)
    if not is_valid:
        raise ValueError(err)
    return cleaned, []


def parse_fasta_content(fasta_text: str) -> Tuple[str, List[AnnotationFeature]]:
    """Parse FASTA-formatted text and return the first sequence."""
    handle = io.StringIO(fasta_text)
    records = list(SeqIO.parse(handle, "fasta"))
    if not records:
        raise ValueError("No FASTA records found in the input.")
    record = records[0]
    seq = clean_sequence(str(record.seq))
    return seq, []


def _extract_annotations(record: SeqRecord, filter_types: Optional[List[str]] = None) -> List[AnnotationFeature]:
    """Extract GenBank feature annotations as AnnotationFeature objects."""
    EXCLUDED_TYPES = {"source"}
    features: List[AnnotationFeature] = []
    
    exons = []
    
    for feat in record.features:
        if feat.type in EXCLUDED_TYPES:
            continue
            
        if feat.type == "exon":
            exons.append(feat)
            
        if filter_types is not None:
            # Check explicit matches or broad UTR matching
            matches = feat.type in filter_types or ("UTR" in filter_types and "UTR" in feat.type)
            if not matches:
                continue
        # Build a label from qualifiers
        label_parts = []
        for qualifier in ("gene", "product", "note", "locus_tag"):
            if qualifier in feat.qualifiers:
                label_parts.append(feat.qualifiers[qualifier][0])
                break
        label = label_parts[0] if label_parts else feat.type

        start = int(feat.location.start)
        end = int(feat.location.end)
        strand = feat.location.strand if feat.location.strand is not None else 1

        features.append(AnnotationFeature(
            type=feat.type,
            label=label,
            start=start,
            end=end,
            strand=strand,
        ))
        
    has_explicit_introns = any(f.type == "intron" for f in features)
    if not has_explicit_introns and exons:
        exons_sorted = sorted(exons, key=lambda f: int(f.location.start))
        exons_by_gene = {}
        for ex in exons_sorted:
            gene = ex.qualifiers.get("gene", ["unknown_gene"])[0]
            exons_by_gene.setdefault(gene, []).append(ex)
            
        for gene, gene_exons in exons_by_gene.items():
            strand = gene_exons[0].location.strand
            for i in range(len(gene_exons) - 1):
                intron_start = int(gene_exons[i].location.end)
                intron_end = int(gene_exons[i+1].location.start)
                
                if intron_start < intron_end:
                    intron_num = len(gene_exons) - 1 - i if strand == -1 else i + 1
                    label = f"{gene} intron {intron_num}"
                    
                    if filter_types is None or "intron" in filter_types:
                        features.append(AnnotationFeature(
                            type="intron",
                            label=label,
                            start=intron_start,
                            end=intron_end,
                            strand=strand if strand is not None else 1
                        ))

    return features


def search_ncbi(query: str, organism: str = "") -> dict:
    """
    Search NCBI for transcripts. Handles NC_ and NM_ detection.
    """
    query = query.strip()
    upper_query = query.upper()
    
    if upper_query.startswith("NC_"):
        return {"type": "chromosome", "results": [{"id": upper_query, "title": f"Chromosome {upper_query}", "organism": organism or "Unknown"}]}
        
    if upper_query.startswith("NM_"):
        return {"type": "transcript", "results": [{"id": query, "title": f"Transcript {query}", "organism": organism or "Unknown"}]}
        
    # Gene name search
    term = f"{query}[Gene Name] AND biomol_mrna[PROP]"
    if organism:
        term += f" AND {organism}[Organism]"
        
    search_handle = Entrez.esearch(db="nucleotide", term=term, retmax=20)
    search_result = Entrez.read(search_handle)
    search_handle.close()
    
    ids = search_result.get("IdList", [])
    if not ids:
        return {"type": "gene", "results": []}
        
    summary_handle = Entrez.esummary(db="nucleotide", id=",".join(ids))
    summaries = Entrez.read(summary_handle)
    summary_handle.close()
    
    results = []
    for s in summaries:
        acc = s.get("Caption", "")
        if not acc:
            continue
        results.append({
            "id": acc,
            "title": s.get("Title", ""),
            "organism": organism or "Unknown"
        })
        
    return {"type": "gene", "results": results}


def fetch_ncbi(accession: str, filter_types: Optional[List[str]] = None) -> Tuple[str, List[AnnotationFeature]]:
    """
    Fetch a GenBank record by accession number.
    Returns (sequence_string, annotation_features).
    """
    accession = accession.strip()

    fetch_handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="gb",
        retmode="text",
    )
    try:
        record = SeqIO.read(fetch_handle, "genbank")
        features = _extract_annotations(record, filter_types)
        
        # Synthesize a promoter feature so it appears in the frontend annotations list.
        # For transcript accessions, the local record often starts at the TSS, so this
        # can be 0-0. The analysis endpoint resolves promoter sequence separately.
        try:
            _, p_start, p_end, is_rev, tss = compute_promoter_region(record, 1000)
            features.append(AnnotationFeature(
                type="promoter",
                label=f"Promoter upstream of TSS",
                start=p_start,
                end=p_end,
                strand=-1 if is_rev else 1
            ))
        except Exception:
            pass # Skip if no gene/mRNA/CDS found to derive a TSS from
            
        seq = clean_sequence(str(record.seq))
        fetch_handle.close()
    except Exception as e:
        fetch_handle.close()
        raise e

    return seq, features


def get_genomic_promoter(accession: str, promoter_length: int = 1000) -> Tuple[str, int, int, bool, int]:
    """
    Attempt to resolve a transcript accession (NM_) to its genomic location (NC_)
    and fetch the upstream promoter region directly.
    
    Returns: (promoter_seq, start, end, is_reverse, tss)
    """
    # 1. Get Gene ID from accession
    search_handle = Entrez.esearch(db="nucleotide", term=accession)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    if not search_results["IdList"]:
        raise ValueError(f"Could not find NCBI record for {accession}")
    
    uid = search_results["IdList"][0]
    
    link_handle = Entrez.elink(dbfrom="nucleotide", db="gene", id=uid)
    link_results = Entrez.read(link_handle)
    link_handle.close()
    
    gene_ids = []
    for link_set in link_results:
        for link in link_set.get("LinkSetDb", []):
            if link["DbTo"] == "gene":
                gene_ids.extend([l["Id"] for l in link["Link"]])
                
    if not gene_ids:
        raise ValueError(f"No associated Gene ID found for {accession}")
        
    # 2. Get Genomic Coords from Gene ID
    gene_id = gene_ids[0]
    summary_handle = Entrez.esummary(db="gene", id=gene_id)
    summary_results = Entrez.read(summary_handle)
    summary_handle.close()
    
    doc_sum = summary_results["DocumentSummarySet"]["DocumentSummary"][0]
    genomic_info = doc_sum.get("GenomicInfo", [])
    if not genomic_info:
         raise ValueError(f"No genomic location information found for Gene ID {gene_id}")
         
    info = genomic_info[0]
    chr_acc = info.get("ChrAccVer")
    chr_start = int(info.get("ChrStart"))
    chr_stop = int(info.get("ChrStop"))
    
    if not chr_acc:
        raise ValueError(f"Missing chromosome accession for Gene ID {gene_id}")
        
    is_reverse = (chr_start > chr_stop)
    tss = chr_start
    
    if is_reverse:
        p_start = tss + 1
        p_end = tss + promoter_length
    else:
        p_start = max(0, tss - promoter_length)
        p_end = max(0, tss - 1)
        
    fetch_start = p_start + 1
    fetch_end = p_end + 1
    
    fetch_handle = Entrez.efetch(
        db="nucleotide",
        id=chr_acc,
        seq_start=fetch_start,
        seq_stop=fetch_end,
        rettype="fasta",
        retmode="text"
    )
    promoter_record = SeqIO.read(fetch_handle, "fasta")
    fetch_handle.close()
    
    promoter_seq = str(promoter_record.seq)
    if is_reverse:
        from Bio.Seq import Seq
        promoter_seq = str(Seq(promoter_seq).reverse_complement())
        
    return clean_sequence(promoter_seq), p_start, p_end, is_reverse, tss


def compute_promoter_region(record, promoter_length: int = 1000) -> tuple:
    """
    Compute the promoter region from a fetched GenBank record.

    Note: Promoter is defined as upstream region relative to TSS and is an approximation.

    Strategy:
      - Finds the primary gene/mRNA/CDS feature to determine TSS and strand.
      - Forward strand (+): TSS = feature.start → promoter = [TSS-length, TSS]
      - Reverse strand (−): TSS = feature.end  → promoter = [TSS, TSS+length]
      - Coordinates are clamped to [0, seq_len].

    Returns:
        (promoter_seq, prom_start, prom_end, is_reverse_strand, tss_position)
    """
    seq = str(record.seq)
    seq_len = len(seq)

    # Find primary feature: gene > mRNA > CDS (in that priority)
    target_feature = None
    for feat_type in ("gene", "mRNA", "CDS"):
        for feat in record.features:
            if feat.type == feat_type:
                target_feature = feat
                break
        if target_feature is not None:
            break

    if target_feature is None:
        raise ValueError(
            "No gene, mRNA, or CDS feature found in the NCBI record. "
            "Cannot determine TSS for promoter extraction."
        )

    strand = target_feature.location.strand
    is_reverse = (strand == -1)

    if is_reverse:
        # Reverse strand: TSS = highest coordinate (end of feature in fwd notation)
        tss = int(target_feature.location.end)
        prom_start = tss
        prom_end = min(seq_len, tss + promoter_length)
    else:
        # Forward strand: TSS = lowest coordinate (start of feature)
        tss = int(target_feature.location.start)
        prom_start = max(0, tss - promoter_length)
        prom_end = tss

    promoter_seq = clean_sequence(seq[prom_start:prom_end])
    return promoter_seq, prom_start, prom_end, is_reverse, tss
