/**
 * TriplexLab TFO Finder - App Logic
 */

// Automatically use the current host for API requests (works locally and on Render)
const API_BASE = window.location.origin;

// State
let state = {
    sequence: "",
    targetSequence: "",
    targetOffset: 0,
    length: 0,
    annotations: [],
    accession: null,
    isNcbi: false,
    selectedRegion: null, // {start, end}
    candidates: [],
    ttsRegions: [],
    isChromosome: false
};

// -------------------------------------------------------------------------
// DOM Elements
// -------------------------------------------------------------------------
const DOM = {
    // Steps
    steps: document.querySelectorAll('.step'),
    panels: document.querySelectorAll('.panel'),
    btnNext1: document.getElementById('btn-next-1'),
    btnBack2: document.getElementById('btn-back-2'),
    btnNext2: document.getElementById('btn-next-2'),
    btnBack3: document.getElementById('btn-back-3'),
    btnGenerate: document.getElementById('btn-generate'),
    btnBack4: document.getElementById('btn-back-4'),
    btnNew: document.getElementById('btn-new-analysis'),

    // Step 1: Input
    tabBtns: document.querySelectorAll('.tab-btn'),
    tabContents: document.querySelectorAll('.tab-content'),
    manualSeq: document.getElementById('manual-seq'),
    manualLength: document.getElementById('manual-length'),
    btnClearManual: document.getElementById('btn-clear-manual'),
    fastaInput: document.getElementById('fasta-file'),
    fastaDrop: document.getElementById('fasta-drop-area'),
    fastaLabel: document.getElementById('fasta-file-name'),
    ncbiId: document.getElementById('ncbi-id'),
    btnFetchNcbi: document.getElementById('btn-fetch-ncbi'),
    ncbiLoading: document.getElementById('ncbi-loading'),
    ncbiOrganism: document.getElementById('ncbi-organism'),
    inputTypeFeedback: document.getElementById('input-type-feedback'),

    // Step 2: Region
    regionRadios: document.querySelectorAll('input[name="region_type"]'),
    customCoords: document.getElementById('custom-coords-inputs'),
    manualStart: document.getElementById('manual-start'),
    manualEnd: document.getElementById('manual-end'),
    cardAnnotations: document.getElementById('card-annotations'),
    annotationsList: document.getElementById('annotations-list'),
    annotationSelect: document.getElementById('annotation-select'),
    annotationFilterCbs: document.querySelectorAll('.ann-filter-cb'),
    promoterInputs: document.getElementById('promoter-inputs'),
    promoterLength: document.getElementById('promoter-length'),
    previewLength: document.getElementById('preview-length'),
    seqPreviewText: document.getElementById('seq-preview-text'),

    // Step 3: Settings
    ttsMinLen: document.getElementById('tts-min-len'),
    valTtsLen: document.getElementById('val-tts-len'),
    ttsPurineRatio: document.getElementById('tts-purine-ratio'),
    valPurineRatio: document.getElementById('val-purine-ratio'),
    lengthCats: document.querySelectorAll('.length-cat'),
    filterUnique: document.getElementById('filter-unique-only'),
    filterMaxTargetMismatches: document.getElementById('filter-max-target-mismatches'),
    filterMaxTfoPurines: document.getElementById('filter-max-tfo-purines'),
    apiError: document.getElementById('api-error'),
    apiErrorMsg: document.getElementById('api-error-msg'),
    genSpinner: document.getElementById('generate-spinner'),

    // Step 4: Results
    resTtsCount: document.getElementById('res-tts-count'),
    resTfoCount: document.getElementById('res-tfo-count'),
    resUniqueCount: document.getElementById('res-unique-count'),
    masterSeqViewer: document.getElementById('master-sequence-viewer'),
    candidatesTbody: document.getElementById('candidates-tbody'),
    noResults: document.getElementById('no-results'),
    btnCopySeq: document.getElementById('btn-copy-seq'),
    btnExportCsv: document.getElementById('btn-export-csv'),
    sortHeaders: document.querySelectorAll('th.sortable')
};

// -------------------------------------------------------------------------
// Navigation & Toast
// -------------------------------------------------------------------------
function goToStep(stepNum) {
    DOM.steps.forEach((s, idx) => {
        if (idx < stepNum) {
            s.classList.add('completed');
            s.classList.remove('active');
        } else if (idx === stepNum - 1) {
            s.classList.add('active');
            s.classList.remove('completed');
        } else {
            s.classList.remove('active', 'completed');
        }
    });

    DOM.panels.forEach((p, idx) => {
        if (idx === stepNum - 1) {
            p.classList.add('active');
        } else {
            p.classList.remove('active');
        }
    });
}

function showToast(msg, isError = false) {
    const toast = document.getElementById('toast');
    toast.textContent = msg;
    toast.style.borderLeftColor = isError ? 'var(--accent-rose)' : 'var(--accent-cyan)';
    toast.classList.add('show');
    setTimeout(() => toast.classList.remove('show'), 3000);
}

// -------------------------------------------------------------------------
// Step 1: Input Handlers
// -------------------------------------------------------------------------
DOM.tabBtns.forEach(btn => {
    btn.addEventListener('click', () => {
        DOM.tabBtns.forEach(b => b.classList.remove('active'));
        DOM.tabContents.forEach(c => c.classList.remove('active'));
        btn.classList.add('active');
        document.getElementById(btn.dataset.tab).classList.add('active');
    });
});

// Manual
DOM.manualSeq.addEventListener('input', (e) => {
    let raw = e.target.value.toUpperCase().replace(/[^ACTGN\s]/g, '');
    let clean = raw.replace(/\s/g, '');
    DOM.manualLength.textContent = `Length: ${clean.length} bp`;
});

DOM.btnClearManual.addEventListener('click', () => {
    DOM.manualSeq.value = '';
    DOM.manualLength.textContent = 'Length: 0 bp';
});

// FASTA Upload
DOM.fastaDrop.addEventListener('dragover', e => { e.preventDefault(); DOM.fastaDrop.classList.add('dragover'); });
DOM.fastaDrop.addEventListener('dragleave', () => DOM.fastaDrop.classList.remove('dragover'));
DOM.fastaDrop.addEventListener('drop', e => {
    e.preventDefault();
    DOM.fastaDrop.classList.remove('dragover');
    if (e.dataTransfer.files.length) {
        DOM.fastaInput.files = e.dataTransfer.files;
        handleFastaUpload();
    }
});
DOM.fastaInput.addEventListener('change', handleFastaUpload);

function handleFastaUpload() {
    const file = DOM.fastaInput.files[0];
    if (file) DOM.fastaLabel.textContent = `Selected: ${file.name}`;
}

// NCBI Fetch & Search
DOM.btnFetchNcbi.addEventListener('click', async () => {
    const id = DOM.ncbiId.value.trim();
    const organism = DOM.ncbiOrganism.value.trim();
    
    if (!id) return showToast('Enter an accession ID or gene name', true);
    
    DOM.btnFetchNcbi.disabled = true;
    DOM.ncbiLoading.classList.remove('hidden');
    document.getElementById('ncbi-loading-text').textContent = 'Searching NCBI...';
    DOM.inputTypeFeedback.textContent = '';
    
    try {
        const res = await fetch(`${API_BASE}/sequence/search?query=${encodeURIComponent(id)}&organism=${encodeURIComponent(organism)}`);
        const data = await res.json();
        
        let accessionToFetch = id;
        
        if (data.type === "error") {
            DOM.inputTypeFeedback.textContent = data.message;
            DOM.inputTypeFeedback.style.color = "var(--accent-rose)";
            return;
        }
        
        if (data.type === "transcript" || data.type === "gene" || data.type === "chromosome") {
            if (!data.results || data.results.length === 0) {
                DOM.inputTypeFeedback.textContent = "No valid records found in NCBI.";
                DOM.inputTypeFeedback.style.color = "var(--accent-rose)";
                return;
            }
            accessionToFetch = data.results[0].id;
        }

        document.getElementById('ncbi-loading-text').textContent = 'Fetching sequence...';
        
        const fetchRes = await fetch(`${API_BASE}/sequence/fetch/${encodeURIComponent(accessionToFetch)}`);
        if (!fetchRes.ok) throw new Error("Could not fetch sequence from NCBI");
        const fetchPayload = await fetchRes.json();
        
        state.sequence = fetchPayload.sequence;
        state.length = fetchPayload.length;
        state.annotations = fetchPayload.annotations || [];
        state.accession = accessionToFetch;
        state.isNcbi = true;
        
        DOM.inputTypeFeedback.textContent = `Auto-Selected: ${accessionToFetch}`;
        DOM.inputTypeFeedback.style.color = "var(--accent-cyan)";
        
        showToast(`Successfully fetched sequence (${fetchPayload.length} bp)`);
        populateRegionStep();
        goToStep(2);
        
    } catch (err) {
        showToast("Error scanning NCBI entries.", true);
    } finally {
        DOM.btnFetchNcbi.disabled = false;
        DOM.ncbiLoading.classList.add('hidden');
    }
});

// Proceed from Step 1 manually or via fasta
DOM.btnNext1.addEventListener('click', async () => {
    const activeTab = document.querySelector('.tab-btn.active').dataset.tab;
    
    if (activeTab === 'tab-manual') {
        const seq = DOM.manualSeq.value.replace(/\s/g, '').toUpperCase();
        if (!seq) return showToast('Please enter a sequence', true);
        
        try {
            const res = await fetch(`${API_BASE}/sequence/parse`, {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({input_type: "manual", sequence: seq})
            });
            if (!res.ok) throw new Error(await res.text());
            const data = await res.json();
            
            state.sequence = data.sequence;
            state.length = data.length;
            state.annotations = [];
            state.isNcbi = false;
            state.accession = null;
            state.isChromosome = false;
            populateRegionStep();
            goToStep(2);
        } catch (e) {
            showToast('Parsing error. Check sequence.', true);
        }
    } 
    else if (activeTab === 'tab-fasta') {
        const file = DOM.fastaInput.files[0];
        if (!file) return showToast('Please select a FASTA file', true);
        
        const reader = new FileReader();
        reader.onload = async (e) => {
            try {
                const res = await fetch(`${API_BASE}/sequence/parse`, {
                    method: 'POST',
                    headers: {'Content-Type': 'application/json'},
                    body: JSON.stringify({input_type: "fasta", fasta_content: e.target.result})
                });
                if (!res.ok) throw new Error(await res.text());
                const data = await res.json();
                
                state.sequence = data.sequence;
                state.length = data.length;
                state.annotations = [];
                state.isNcbi = false;
                state.accession = null;
                state.isChromosome = false;
                populateRegionStep();
                goToStep(2);
            } catch (err) {
                showToast('Invalid FASTA format', true);
            }
        };
        reader.readAsText(file);
    }
});


// -------------------------------------------------------------------------
// Step 2: Region Selection Handlers
// -------------------------------------------------------------------------
function populateRegionStep() {
    DOM.previewLength.textContent = state.length;
    DOM.seqPreviewText.textContent = state.length > 500 
        ? state.sequence.substring(0, 500) + '... [truncated]' 
        : state.sequence;
    
    DOM.manualStart.max = state.length;
    DOM.manualEnd.max = state.length;
    DOM.manualEnd.value = state.length;
    
    // Annotations UI
    if (state.annotations.length > 0) {
        DOM.cardAnnotations.style.display = 'block';
        updateAnnotationsDropdown();
    } else {
        DOM.cardAnnotations.style.display = 'none';
        if (document.getElementById('region-annotation').checked) {
            document.getElementById('region-manual').checked = true;
        }
    }

    const fullRegionCard = document.getElementById('region-full').closest('.radio-card');
    if (state.isNcbi) {
        document.getElementById('region-full').disabled = true;
        fullRegionCard.style.opacity = '0.4';
        fullRegionCard.style.pointerEvents = 'none';

        if (state.annotations.length > 0) {
            document.getElementById('region-annotation').checked = true;
        } else {
            document.getElementById('region-manual').checked = true;
        }
    } else {
        document.getElementById('region-full').disabled = false;
        fullRegionCard.style.opacity = '1';
        fullRegionCard.style.pointerEvents = 'auto';
        if (!document.querySelector('input[name="region_type"]:checked') || document.getElementById('region-annotation').checked && state.annotations.length === 0) {
            document.getElementById('region-full').checked = true;
        }
    }
    
    updateRegionUI();
}

DOM.regionRadios.forEach(r => {
    r.addEventListener('change', updateRegionUI);
});

function getActiveAnnotationFilters() {
    return Array.from(DOM.annotationFilterCbs)
        .filter(cb => cb.checked)
        .map(cb => cb.value);
}

function updateAnnotationsDropdown() {
    DOM.annotationSelect.innerHTML = '';
    const activeFilters = getActiveAnnotationFilters();
    
    // Create a mapping to easily find the true index later
    state.filteredAnnotations = [];
    
    state.annotations.forEach((ann, originalIndex) => {
        let matches = false;
        if (activeFilters.includes(ann.type)) matches = true;
        if (ann.type.includes("UTR") && activeFilters.includes("UTR")) matches = true;
        
        if (matches) {
            state.filteredAnnotations.push({ ...ann, originalIndex });
            const opt = document.createElement('option');
            opt.value = state.filteredAnnotations.length - 1; // Index in the filtered array
            if (ann.type === 'promoter') {
                opt.textContent = `${ann.label} (computed during analysis) ${ann.strand === -1 ? '(rev)' : '(fwd)'}`;
            } else {
                opt.textContent = `${ann.type}: ${ann.label} (${ann.start}-${ann.end}) ${ann.strand === -1 ? '(rev)' : '(fwd)'}`;
            }
            DOM.annotationSelect.appendChild(opt);
        }
    });
    
    if (state.filteredAnnotations.length === 0) {
        const opt = document.createElement('option');
        opt.value = "";
        opt.textContent = "No annotations match filters";
        opt.disabled = true;
        opt.selected = true;
        DOM.annotationSelect.appendChild(opt);
    }

    updatePromoterLengthVisibility();
}

DOM.annotationFilterCbs.forEach(cb => {
    cb.addEventListener('change', () => {
        updateAnnotationsDropdown();
        updateRegionUI();
    });
});

DOM.annotationSelect.addEventListener('change', () => {
    updatePromoterLengthVisibility();
});

function updatePromoterLengthVisibility() {
    const selIdx = parseInt(DOM.annotationSelect.value);
    const isPromoterSelected = !isNaN(selIdx) &&
        state.filteredAnnotations &&
        state.filteredAnnotations[selIdx] &&
        state.filteredAnnotations[selIdx].type === 'promoter';

    if (isPromoterSelected) {
        DOM.promoterInputs.classList.remove('hidden');
    } else {
        DOM.promoterInputs.classList.add('hidden');
    }
}

function updateRegionUI() {
    document.querySelectorAll('.radio-card').forEach(c => c.classList.remove('selected'));
    const checked = document.querySelector('input[name="region_type"]:checked');
    if (checked) checked.closest('.radio-card').classList.add('selected');
    
    DOM.customCoords.classList.add('hidden');
    DOM.annotationsList.classList.add('hidden');
    DOM.promoterInputs.classList.add('hidden');

    if (checked && checked.value === 'manual') DOM.customCoords.classList.remove('hidden');
    if (checked && checked.value === 'annotation') {
        DOM.annotationsList.classList.remove('hidden');
        updatePromoterLengthVisibility();
    }
}

DOM.btnBack2.addEventListener('click', () => goToStep(1));
DOM.btnNext2.addEventListener('click', () => {
    const type = document.querySelector('input[name="region_type"]:checked').value;
    
    if (type === 'full') {
        state.selectedRegion = null;
        state.regionType = 'full';
    } else if (type === 'manual') {
        const s = parseInt(DOM.manualStart.value) || 1;
        const e = parseInt(DOM.manualEnd.value) || state.length;
        if (s >= e || s < 1 || e > state.length) {
            return showToast('Invalid custom coordinates', true);
        }
        state.selectedRegion = { start: s - 1, end: e }; // 0-indexed for backend
        state.regionType = 'manual';
    } else if (type === 'annotation') {
        const selIdx = parseInt(DOM.annotationSelect.value);
        if (isNaN(selIdx) || !state.filteredAnnotations || !state.filteredAnnotations[selIdx]) {
            return showToast('Please select a valid annotation', true);
        }
        const ann = state.filteredAnnotations[selIdx];
        if (ann.type === 'promoter') {
            state.selectedRegion = null;
            state.regionType = 'promoter';
            state.promoterLength = parseInt(DOM.promoterLength.value) || 1000;
        } else {
            state.selectedRegion = { start: ann.start, end: ann.end, label: ann.label, feature_type: ann.type };
            state.regionType = 'annotation';
        }
    }
    
    goToStep(3);
});


// -------------------------------------------------------------------------
// Step 3: Settings & API Call
// -------------------------------------------------------------------------
DOM.ttsMinLen.addEventListener('input', e => DOM.valTtsLen.textContent = e.target.value + ' bp');
DOM.ttsPurineRatio.addEventListener('input', e => DOM.valPurineRatio.textContent = Math.round(e.target.value * 100) + '%');

DOM.btnBack3.addEventListener('click', () => goToStep(2));
DOM.btnGenerate.addEventListener('click', async () => {
    DOM.apiError.classList.add('hidden');
    DOM.btnGenerate.disabled = true;
    DOM.genSpinner.classList.remove('hidden');
    document.querySelector('#btn-generate .btn-text').textContent = 'Processing...';
    
    // Gather length categories
    const lengths = Array.from(DOM.lengthCats)
        .filter(cb => cb.checked)
        .map(cb => cb.value);
        
    if (lengths.length === 0) {
        DOM.apiError.classList.remove('hidden');
        DOM.apiErrorMsg.textContent = "Please select at least one TFO length category.";
        resetGenerateBtn();
        return;
    }

    const payload = {
        sequence: state.isNcbi ? "" : state.sequence,
        accession: state.isNcbi ? state.accession : undefined,
        region: state.selectedRegion,
        region_type: state.regionType || 'full',
        promoter_length: state.promoterLength || 1000,
        is_ncbi: state.isNcbi || false,
        tts_min_length: parseInt(DOM.ttsMinLen.value),
        tts_purine_ratio: parseFloat(DOM.ttsPurineRatio.value),
        filters: {
            only_unique: DOM.filterUnique.checked,
            max_target_mismatches: parseInt(DOM.filterMaxTargetMismatches.value) || 0,
            max_tfo_purines: parseInt(DOM.filterMaxTfoPurines.value) >= 0 ? parseInt(DOM.filterMaxTfoPurines.value) : 100,
            length_categories: lengths
        }
    };

    try {
        const res = await fetch(`${API_BASE}/tfo/find`, {
            method: 'POST',
            headers: {'Content-Type': 'application/json'},
            body: JSON.stringify(payload)
        });
        
        if (!res.ok) throw new Error(await res.text());
        const data = await res.json();
        
        state.candidates = data.candidates || [];
        state.ttsRegions = data.tts_regions || [];
        state.targetSequence = data.target_sequence || state.sequence || "";
        state.targetOffset = data.target_offset || 0;
        
        renderResults();
        goToStep(4);
    } catch (err) {
        DOM.apiError.classList.remove('hidden');
        let errorMsg = "An error occurred during pipeline execution.";
        try {
            const parsed = JSON.parse(err.message);
            errorMsg = parsed.detail || errorMsg;
        } catch (e) {
            errorMsg = err.message || errorMsg;
        }
        DOM.apiErrorMsg.textContent = errorMsg;
    } finally {
        resetGenerateBtn();
    }
});

function resetGenerateBtn() {
    DOM.btnGenerate.disabled = false;
    DOM.genSpinner.classList.add('hidden');
    document.querySelector('#btn-generate .btn-text').textContent = 'Find TFOs';
}


// -------------------------------------------------------------------------
// Step 4: Results Render & Sorting
// -------------------------------------------------------------------------
let sortCol = 'score';
let sortDesc = true;

function renderResults() {
    // Stats
    DOM.resTtsCount.textContent = state.ttsRegions.length;
    DOM.resTfoCount.textContent = state.candidates.length;
    DOM.resUniqueCount.textContent = state.candidates.filter(c => c.is_unique).length;
    DOM.btnExportCsv.disabled = state.candidates.length === 0;

    // Build Master Sequence Viewer
    let seqHTML = "";
    let lastIdx = 0;
    const fullSeq = state.targetSequence || state.sequence || "";
    const targetOffset = state.targetOffset || 0;
    
    // Sort TTS regions purely for display (sequential)
    const sortedTTS = [...state.ttsRegions].sort((a,b) => a.start - b.start);
    
    // To prevent browser lag on huge sequences, we'll only render 
    // the regions containing highlights plus some padding.
    if (fullSeq.length > 50000 && sortedTTS.length > 0) {
        seqHTML = `<i>Sequence too large to display in full (${(fullSeq.length/1e6).toFixed(1)} MB). Showing target regions:</i><br><br>`;
        sortedTTS.forEach(tts => {
            const localStart = Math.max(0, tts.start - targetOffset);
            const localEnd = Math.min(fullSeq.length - 1, tts.end - targetOffset);
            const contextStart = Math.max(0, localStart - 50);
            const contextEnd = Math.min(fullSeq.length, localEnd + 51);
            seqHTML += `<div class="seq-context">...${fullSeq.substring(contextStart, localStart)}<span class="tts-highlight" title="TTS [${tts.start}-${tts.end}]">${fullSeq.substring(localStart, localEnd + 1)}</span>${fullSeq.substring(localEnd + 1, contextEnd)}...</div>`;
        });
    } else {
        sortedTTS.forEach(tts => {
            const localStart = Math.max(0, tts.start - targetOffset);
            const localEnd = Math.min(fullSeq.length - 1, tts.end - targetOffset);
            if (localStart > lastIdx) {
                seqHTML += fullSeq.substring(lastIdx, localStart);
            }
            seqHTML += `<span class="tts-highlight" title="TTS [${tts.start}-${tts.end}] Length: ${tts.length}">${fullSeq.substring(localStart, localEnd + 1)}</span>`;
            lastIdx = localEnd + 1;
        });
        if (lastIdx < fullSeq.length) {
            seqHTML += fullSeq.substring(lastIdx);
        }
    }
    DOM.masterSeqViewer.innerHTML = seqHTML || "No sequence loaded.";

    renderTable();
}

function renderTable() {
    DOM.candidatesTbody.innerHTML = '';
    
    if (state.candidates.length === 0) {
        DOM.noResults.classList.remove('hidden');
        return;
    }
    DOM.noResults.classList.add('hidden');
    
    // Sort natively
    const sorted = [...state.candidates].sort((a, b) => {
        let valA = a[sortCol];
        let valB = b[sortCol];
        if (sortDesc) return valB > valA ? 1 : -1;
        return valA > valB ? 1 : -1;
    });

    sorted.forEach((c, i) => {
        const tr = document.createElement('tr');
        
        // Match color logic
        let matchClass = 'color-match';
        if (c.match_count > 1) matchClass += ' warn';
        
        const purines = c.purine_count ?? 0;
        
        tr.innerHTML = `
            <td style="color:var(--text-secondary);font-size:0.8rem">${c.start}-${c.end}</td>
            <td style="letter-spacing:1px;font-weight:500;">${c.tfo_sequence || ''}</td>
            <td>${c.length || 0}</td>
            <td>${c.target_mismatches ?? 0}</td>
            <td style="font-weight:bold;color:var(--accent-cyan)">${purines}</td>
            <td class="${matchClass}">${c.match_count ?? 0} ${c.is_unique ? 'Unique' : ''}</td>
        `;
        DOM.candidatesTbody.appendChild(tr);
    });
}

// Table Sorting
DOM.sortHeaders.forEach(th => {
    th.addEventListener('click', () => {
        const col = th.dataset.sort;
        if (sortCol === col) {
            sortDesc = !sortDesc; // toggle
        } else {
            sortCol = col;
            sortDesc = (col === 'length' || col === 'tfo_sequence' || col === 'score'); // Desc for these, Asc for matches/mismatches
        }
        
        // Update header UI
        DOM.sortHeaders.forEach(h => h.classList.remove('desc', 'asc'));
        th.classList.add(sortDesc ? 'desc' : 'asc');
        
        renderTable();
    });
});

// CSV Export
DOM.btnExportCsv.addEventListener('click', () => {
    if(!state.candidates.length) return;
    const headers = ["Start", "End", "Length", "Category", "TTS_Sequence", "TFO_Sequence", "Target_Mismatches", "Purines", "Binding_Occurrences", "Is_Unique", "Score"];
    const rows = state.candidates.map(c => [
        c.start, c.end, c.length, c.length_category,
        c.tts_sequence, c.tfo_sequence, c.target_mismatches, c.purine_count,
        c.match_count, c.is_unique, c.score
    ]);
    
    let csv = headers.join(',') + '\n';
    rows.forEach(r => { csv += r.join(',') + '\n'; });
    
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'TFO_Candidates.csv';
    a.click();
    window.URL.revokeObjectURL(url);
});

DOM.btnCopySeq.addEventListener('click', () => {
    navigator.clipboard.writeText(state.sequence);
    showToast('Full target sequence copied');
});

DOM.btnBack4.addEventListener('click', () => goToStep(3));
DOM.btnNew.addEventListener('click', () => {
    state = { sequence: "", targetSequence: "", targetOffset: 0, length: 0, annotations: [], accession: null, isNcbi: false, selectedRegion: null, regionType: 'full', promoterLength: 1000, candidates: [], ttsRegions: [], isChromosome: false };
    DOM.manualSeq.value = '';
    DOM.manualLength.textContent = 'Length: 0 bp';
    DOM.fastaInput.value = '';
    DOM.fastaLabel.textContent = '';
    DOM.ncbiId.value = '';
    goToStep(1);
});
