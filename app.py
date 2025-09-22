import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import io
import numpy as np
from collections import Counter
from Bio import Entrez, SeqIO

from utils import (
    parse_fasta, gc_content, reverse_complement, wrap
)
from motifs import get_basic_stats
from motifs.base_motif import select_best_nonoverlapping_motifs
# Import the new unified orchestrator
from all_motifs_refactored import all_motifs_refactored
from motifs import visualization as viz
from motifs.enhanced_visualization import create_comprehensive_information_based_visualizations
# Import class definitions for selection interface
from class_definitions import CLASS_DEFINITIONS, DEFAULT_SELECTED_CLASSES, DEFAULT_SELECTED_SUBCLASSES

# Import new configuration modules
try:
    from classification_config import MOTIF_LENGTH_LIMITS, SCORING_METHODS, get_motif_limits
    CONFIG_AVAILABLE = True
except ImportError:
    CONFIG_AVAILABLE = False
    MOTIF_LENGTH_LIMITS = {}
    SCORING_METHODS = {}

# ---------- PAGE CONFIG ----------
st.set_page_config(
    page_title="Non-B DNA Motif Finder",
    layout="wide",
    page_icon="🧬",
    menu_items={'About': "Non-B DNA Motif Finder | Developed by Dr. Venkata Rajesh Yella"}
)


# ---------- PATCH: Ensure every motif has Subclass ----------
def ensure_subclass(motif):
    """Guarantee every motif has a string 'Subclass'"""
    if isinstance(motif, dict):
        if 'Subclass' not in motif or motif['Subclass'] is None:
            motif['Subclass'] = motif.get('Subtype', 'Other')
        return motif
    else:
        # Handle non-dict motifs gracefully (could log/warn here)
        return {'Subclass': 'Other', 'Motif': motif}


# ---------- PROFESSIONAL CSS FOR BALANCED DESIGN ----------
st.markdown("""
    <style>
    body, [data-testid="stAppViewContainer"], .main {
        background: #f7fafd !important;
        font-family: 'Montserrat', Arial, sans-serif !important;
    }
    /* Tabs: medium-large, bold, clean */
    .stTabs [data-baseweb="tab-list"] {
        width: 120vw !important;
        justify-content: stretch !important;
        border-bottom: 2px solid #1565c0;
        background: linear-gradient(90deg,#eaf3fa 0%,#f7fafd 100%) !important;
        box-shadow: 0 2px 8px #dae5f2;
        margin-bottom: 0;
    }
    .stTabs [data-baseweb="tab"] {
        font-size: 1.75rem !important;
        font-weight: 700 !important;
        flex: 1 1 0%;
        min-width: 0 !important;
        padding: 15px 0 15px 0 !important;
        text-align: center;
        color: #1565c0 !important;
        background: #eaf3fa !important;
        border-right: 1px solid #eee !important;
        letter-spacing: 0.03em;
    }
    .stTabs [aria-selected="true"] {
        color: #002147 !important;
        border-bottom: 5px solid #1565c0 !important;
        background: #f7fafd !important;
        box-shadow: 0 4px 8px #e0e5ea;
    }
    .stTabs [data-baseweb="tab"]:last-child {
        border-right: none !important;
    }
    /* Headings: harmonized medium size */
    h1, h2, h3, h4 {
        font-family: 'Montserrat', Arial, sans-serif !important;
        color: #1565c0 !important;
        font-weight: 800 !important;
        margin-top: 0.8em;
        margin-bottom: 0.8em;
    }
    h1 { font-size:2.05rem !important; }
    h2 { font-size:1.55rem !important; }
    h3 { font-size:1.19rem !important; }
    h4 { font-size:1.09rem !important; }
    /* Markdown/text: medium size, easy reading with proper spacing */
    .stMarkdown, .markdown-text-container, .stText, p, span, label, input, .stTextInput>div>div>input, .stSelectbox>div>div>div, .stMultiSelect>div>div>div, .stRadio>div>div>label>div {
        font-size: 1.08rem !important;
        font-family: 'Montserrat', Arial, sans-serif !important;
        line-height: 1.6 !important;
        margin-bottom: 0.5rem !important;
    }
    /* Buttons: modern, medium */
    .stButton>button {
        font-size: 1.08rem !important;
        font-family: 'Montserrat', Arial, sans-serif !important;
        padding: 0.45em 1.2em !important;
        background: linear-gradient(90deg,#1565c0 0%,#2e8bda 100%) !important;
        color: #fff !important;
        border-radius: 7px !important;
        border: none !important;
        font-weight: 600 !important;
        box-shadow: 0 2px 8px #b5cbe6;
        transition: background 0.2s;
    }
    .stButton>button:hover {
        background: linear-gradient(90deg,#2e8bda 0%,#1565c0 100%) !important;
    }
    /* DataFrame font with better spacing */
    .stDataFrame, .stTable {
        font-size: 1.05rem !important;
        font-family: 'Montserrat', Arial, sans-serif !important;
    }
    /* Improved tab content spacing to prevent overlap */
    .stTabs [data-baseweb="tab-panel"] {
        padding-top: 2rem !important;
        padding-left: 1rem !important;
        padding-right: 1rem !important;
    }
    /* Better spacing for analysis summary cards */
    .analysis-summary-card {
        margin: 1.5rem 0 !important;
        padding: 1.5rem !important;
    }
    /* Fix potential text overlap in selectbox and multiselect */
    .stSelectbox > div > div > div, .stMultiSelect > div > div > div {
        min-height: 3.0rem !important;
        padding: 0.8rem !important;
        line-height: 1.5 !important;
    }
    /* Enhanced input field spacing to prevent placeholder overlap */
    .stTextInput > div > div > input, .stTextArea textarea {
        padding: 0.75rem 1rem !important;
        min-height: 2.8rem !important;
        line-height: 1.4 !important;
    }
    /* Better spacing for radio button options */
    .stRadio > div > div > label {
        margin-bottom: 0.5rem !important;
        padding: 0.3rem 0.5rem !important;
    }
    /* Number input spacing */
    .stNumberInput > div > div > input {
        padding: 0.75rem 1rem !important;
        min-height: 2.8rem !important;
    }
    </style>
""", unsafe_allow_html=True)


# ---------- CONSTANTS ----------
# Use classification config if available, otherwise fallback to defaults
if CONFIG_AVAILABLE:
    MOTIF_ORDER = list(MOTIF_LENGTH_LIMITS.keys())
    # Expand special cases
    expanded_order = []
    for motif in MOTIF_ORDER:
        if motif == "Slipped_DNA_DR":
            expanded_order.extend(["Slipped DNA (Direct Repeat)", "Slipped DNA (STR)"])
        elif motif == "Slipped_DNA_STR":
            continue  # Already added above
        elif motif == "eGZ":
            expanded_order.append("eGZ (Extruded-G)")
        elif motif == "G4":
            expanded_order.extend(["G4", "Relaxed G4", "Bulged G4", "Bipartite G4", "Multimeric G4"])
        elif motif == "G-Triplex":
            expanded_order.append("G-Triplex")
        elif motif == "AC-motif":
            expanded_order.append("AC-Motif")
        elif motif == "A-philic_DNA":
            expanded_order.append("A-philic DNA")  # NEW: Class 9
        else:
            # Map to display names
            display_name = motif.replace("_", " ").replace("-", "-")
            if motif == "Curved_DNA":
                display_name = "Curved DNA"
            elif motif == "Z-DNA":
                display_name = "Z-DNA"
            elif motif == "R-Loop":
                display_name = "R-Loop"
            elif motif == "Triplex":
                display_name = "Triplex DNA"
            elif motif == "Sticky_DNA":
                display_name = "Sticky DNA"
            elif motif == "i-Motif":
                display_name = "i-Motif"
            expanded_order.append(display_name)
    
    MOTIF_ORDER = expanded_order + ["Hybrid", "Non-B DNA Clusters"]  # Classes 10, 11
else:
    # Fallback to original order with A-philic DNA added
    MOTIF_ORDER = [
        "Sticky DNA","Curved DNA","Z-DNA","eGZ (Extruded-G)","Slipped DNA","R-Loop",
        "Cruciform","Triplex DNA","G-Triplex","G4","Relaxed G4","Bulged G4","Bipartite G4",
        "Multimeric G4","i-Motif","AC-Motif","A-philic DNA","Hybrid","Non-B DNA Clusters"
    ]

MOTIF_COLORS = {
    "Curved DNA": "#FF9AA2","Z-DNA": "#FFB7B2","eGZ (Extruded-G)": "#6A4C93",
    "Slipped DNA": "#FFDAC1","Slipped DNA (Direct Repeat)": "#FFDAC1","Slipped DNA (STR)": "#FFE4B3",
    "R-Loop": "#FFD3B6","Cruciform": "#E2F0CB",
    "Triplex DNA": "#B5EAD7","Sticky DNA": "#DCB8CB","G-Triplex": "#C7CEEA",
    "G4": "#A2D7D8","Relaxed G4": "#A2D7B8","Bulged G4": "#A2A7D8",
    "Bipartite G4": "#A2D788","Multimeric G4": "#A2A7B8","i-Motif": "#B0C4DE",
    "A-philic DNA": "#E6B8F7",  # NEW: Light purple for A-philic DNA
    "Hybrid": "#C1A192","Non-B DNA Clusters": "#A2C8CC","AC-Motif": "#F5B041"
}
PAGES = {
    "Home": "Overview",
    "Upload & Analyze": "Sequence Upload and Motif Analysis",
    "Results": "Analysis Results and Visualization",
    "Download": "Export Data",
    "Documentation": "Scientific Documentation & References"
}
Entrez.email = "raazbiochem@gmail.com"
Entrez.api_key = None

EXAMPLE_FASTA = """>Example Sequence
ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC
ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT
GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA
GAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA
CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT
"""
EXAMPLE_MULTI_FASTA = """>Seq1
ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC
ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT
>Seq2
GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA
GAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG
>Seq3
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA
"""

# Streamlined session state - removed unnecessary configuration switches
for k, v in {
    'seqs': [],
    'names': [],
    'results': [],
    'summary_df': pd.DataFrame(),
    'hotspots': [],
    'analysis_status': "Ready",
    'selected_classes': DEFAULT_SELECTED_CLASSES.copy(),
    'selected_subclasses': DEFAULT_SELECTED_SUBCLASSES.copy()
}.items():
    if k not in st.session_state:
        st.session_state[k] = v

# ---- TABS ----
tabs = st.tabs(list(PAGES.keys()))
tab_pages = dict(zip(PAGES.keys(), tabs))

with tab_pages["Home"]:
    st.markdown("<h1>Non-B DNA Motif Finder</h1>", unsafe_allow_html=True)
    left, right = st.columns([1,1])
    with left:
        try:
            st.image("nbdcircle.JPG")
        except:
            # If image doesn't exist, show placeholder
            st.markdown("""
            <div style='background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                        border-radius: 15px; padding: 40px; text-align: center; color: white;'>
                <h2 style='margin: 0; color: white;'>🧬</h2>
                <h3 style='margin: 10px 0 0 0; color: white;'>NBD Finder</h3>
                <p style='margin: 5px 0 0 0; color: #E8E8E8;'>Non-B DNA Detection</p>
            </div>
            """, unsafe_allow_html=True)
    with right:
        st.markdown("""
        <div style='font-family:Montserrat, Arial; font-size:1.14rem; color:#222; line-height:1.7; padding:18px; background:#f8f9fa; border-radius:14px; box-shadow:0 2px 8px #eee;'>
        <b>Non-canonical DNA structures</b> play key roles in genome stability, regulation, and evolution.<br>
        This application detects and analyzes <b>11 major classes with 22+ subclasses</b> of Non-B DNA motifs in any DNA sequence or multi-FASTA file.<br>
        <b>Motif Classes (11 classes, 22+ subclasses):</b><br>
        <span style='color:#1565c0;'>
            <b>1. Curved DNA</b> (Global curvature, Local Curvature),<br>
            <b>2. Slipped DNA</b> (Direct Repeat, STR),<br>
            <b>3. Cruciform DNA</b> (Inverted Repeats),<br>
            <b>4. R-loop</b> (R-loop formation sites),<br>
            <b>5. Triplex</b> (Triplex, Sticky DNA),<br>
            <b>6. G-Quadruplex Family</b> (Multimeric G4, Canonical G4, Relaxed G4, Bulged G4, Bipartite G4, Imperfect G4, G-Triplex intermediate),<br>
            <b>7. i-Motif Family</b> (Canonical i-motif, Relaxed i-motif, AC-motif),<br>
            <b>8. Z-DNA</b> (Z-DNA, eGZ (Extruded-G) DNA),<br>
            <b>9. A-philic DNA</b> (A-philic DNA),<br>
            <b>10. Hybrid</b> (dynamic overlaps),<br>
            <b>11. Non-B DNA Clusters</b> (dynamic clusters).
        </span>
        <br>
        <b>Upload single or multi-FASTA files...</b>
        </div>
        """, unsafe_allow_html=True)

# Updated Streamlit layout: Input Method + Sequence Preview + Analysis side-by-side
#  "Upload & Analyze" 
# It expects helper functions and constants to exist in the module:
# parse_fasta, get_basic_stats, EXAMPLE_FASTA, EXAMPLE_MULTI_FASTA, parse_fasta, all_motifs_refactored,
# ensure_subclass, MOTIF_ORDER, wrap, Entrez, SeqIO, Counter, pd, st

with tab_pages["Upload & Analyze"]:
    st.markdown("<h2>Sequence Upload and Motif Analysis</h2>", unsafe_allow_html=True)
    st.markdown('<span style="font-family:Montserrat,Arial; font-size:1.12rem;">Supports multi-FASTA and single FASTA. Paste, upload, select example, or fetch from NCBI.</span>', unsafe_allow_html=True)
    st.caption("Supported formats: .fa, .fasta, .txt | Limit: 200MB/file.")

    st.markdown("---")

    # Create two main columns: left = Input & Preview, right = Analysis controls & Summary
    col_left, col_right = st.columns([1, 1], gap="large")

    # ----- LEFT COLUMN: Input Method + Sequence Preview -----
    with col_left:
        st.markdown("### 📁 Input Method")
        input_method = st.radio("Choose your input method:",
                                ["📂 Upload FASTA File", "✏️ Paste Sequence", "🧪 Example Data", "🌐 NCBI Fetch"],
                                horizontal=True)

        seqs, names = [], []

        if input_method == "📂 Upload FASTA File":
            fasta_file = st.file_uploader("Drag and drop FASTA/multi-FASTA file here", type=["fa", "fasta", "txt"])
            if fasta_file:
                content = fasta_file.read().decode("utf-8")
                seqs, names = [], []
                cur_seq, cur_name = "", ""
                for line in content.splitlines():
                    if line.startswith(">"):
                        if cur_seq:
                            seqs.append(parse_fasta(cur_seq))
                            names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                        cur_name = line.strip().lstrip(">")
                        cur_seq = ""
                    else:
                        cur_seq += line.strip()
                if cur_seq:
                    seqs.append(parse_fasta(cur_seq))
                    names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                if seqs:
                    st.success(f"✅ Loaded {len(seqs)} sequences.")
                    for i, seq in enumerate(seqs[:3]):
                        stats = get_basic_stats(seq)
                        st.markdown(f"**{names[i]}**: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                        st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                    if len(seqs) > 3:
                        st.caption(f"...and {len(seqs)-3} more.")
                else:
                    st.warning("No sequences found.")

        elif input_method == "✏️ Paste Sequence":
            seq_input = st.text_area("Paste single or multi-FASTA here:", height=150, placeholder="Paste your DNA sequence(s) here...")
            if seq_input:
                seqs, names = [], []
                cur_seq, cur_name = "", ""
                for line in seq_input.splitlines():
                    if line.startswith(">"):
                        if cur_seq:
                            seqs.append(parse_fasta(cur_seq))
                            names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                        cur_name = line.strip().lstrip(">")
                        cur_seq = ""
                    else:
                        cur_seq += line.strip()
                if cur_seq:
                    seqs.append(parse_fasta(cur_seq))
                    names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                if seqs:
                    st.success(f"✅ Pasted {len(seqs)} sequences.")
                    for i, seq in enumerate(seqs[:3]):
                        stats = get_basic_stats(seq)
                        st.markdown(f"**{names[i]}**: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                        st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                    if len(seqs) > 3:
                        st.caption(f"...and {len(seqs)-3} more.")
                else:
                    st.warning("No sequences found.")

        elif input_method == "🧪 Example Data":
            ex_type = st.radio("Example Type:", ["Single Example", "Multi-FASTA Example"], horizontal=True)
            if ex_type == "Single Example":
                if st.button("🔬 Load Single Example"):
                    seqs = [parse_fasta(EXAMPLE_FASTA)]
                    names = ["Example Sequence"]
                    st.success("✅ Single example sequence loaded.")
                    stats = get_basic_stats(seqs[0])
                    st.code(EXAMPLE_FASTA, language="fasta")
                    st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
            else:
                if st.button("🔬 Load Multi-FASTA Example"):
                    seqs, names = [], []
                    cur_seq, cur_name = "", ""
                    for line in EXAMPLE_MULTI_FASTA.splitlines():
                        if line.startswith(">"):
                            if cur_seq:
                                seqs.append(parse_fasta(cur_seq))
                                names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                            cur_name = line.strip().lstrip(">")
                            cur_seq = ""
                        else:
                            cur_seq += line.strip()
                    if cur_seq:
                        seqs.append(parse_fasta(cur_seq))
                        names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                    st.success(f"✅ Multi-FASTA example loaded with {len(seqs)} sequences.")
                    for i, seq in enumerate(seqs[:3]):
                        stats = get_basic_stats(seq)
                        st.markdown(f"**{names[i]}**: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                        st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                    st.code(EXAMPLE_MULTI_FASTA, language="fasta")

        elif input_method == "🌐 NCBI Fetch":
            db = st.radio("NCBI Database", ["nucleotide", "gene"], horizontal=True,
                          help="Only nucleotide and gene databases are applicable for DNA motif analysis")
            query_type = st.radio("Query Type", ["Accession", "Gene Name", "Custom Query"], horizontal=True)
            motif_examples = {
                "G-quadruplex": "NR_003287.2 (human telomerase RNA)",
                "Z-DNA": "NM_001126112.2 (human ADAR1 gene)",
                "R-loop": "NR_024540.1 (human SNRPN gene)",
                "eGZ-motif": "CGG repeat region",
                "AC-motif": "A-rich/C-rich consensus region"
            }
            with st.expander("Motif Example Queries"):
                for motif, example in motif_examples.items():
                    st.write(f"**{motif}**: `{example}`")
            query = st.text_input("Enter query (accession, gene, etc.):")
            rettype = st.radio("Return Format", ["fasta", "gb"], horizontal=True)
            retmax = st.number_input("Max Records", min_value=1, max_value=20, value=3)
            if st.button("Fetch from NCBI"):
                if query:
                    with st.spinner("Contacting NCBI..."):
                        handle = Entrez.efetch(db=db, id=query, rettype=rettype, retmode="text")
                        records = list(SeqIO.parse(handle, rettype))
                        handle.close()
                        seqs = [str(rec.seq).upper().replace("U", "T") for rec in records]
                        names = [rec.id for rec in records]
                    if seqs:
                        st.success(f"Fetched {len(seqs)} sequences.")
                        for i, seq in enumerate(seqs[:3]):
                            stats = get_basic_stats(seq)
                            st.markdown(f"<b>{names[i]}</b>: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                            st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                else:
                    st.warning("Enter a query before fetching.")

        # Persist sequences to session state if any found from input
        if seqs:
            st.session_state.seqs = seqs
            st.session_state.names = names
            st.session_state.results = []

        # Sequence Preview lives under the input column for immediate feedback
        if st.session_state.get('seqs'):
            st.markdown("### 📊 Sequence Preview")
            for i, seq in enumerate(st.session_state.seqs[:2]):
                stats = get_basic_stats(seq)
                st.markdown(f"**{st.session_state.names[i]}** ({len(seq):,} bp) | GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}", unsafe_allow_html=True)
                st.code(wrap(seq[:400]), language="fasta")
            if len(st.session_state.seqs) > 2:
                st.caption(f"...and {len(st.session_state.seqs)-2} more.")

    # ----- RIGHT COLUMN: Analysis Controls + Run Button + Summary Table -----
    with col_right:
        st.markdown("### 🚀 Analysis & Run")
        
        # ========== SUCCINCT MOTIF CLASS SELECTOR ========== 
        
        # Quick action buttons
        qa1, qa2, qa3 = st.columns([1,1,1])
        with qa1:
            if st.button("🎯 All", use_container_width=True, help="Select all classes", key="select_all_classes"):
                st.session_state["selected_classes"] = list(CLASS_DEFINITIONS.keys())
        with qa2:
            if st.button("❌ Clear", use_container_width=True, help="Clear selection", key="clear_all_classes"):
                st.session_state["selected_classes"] = []
        with qa3:
            if st.button("🎲 Core", use_container_width=True, help="Select core classes", key="select_core_classes"):
                core = ["Curved_DNA", "G-Quadruplex", "Z-DNA", "Slipped_DNA", "Cruciform"]
                st.session_state["selected_classes"] = [c for c in core if c in CLASS_DEFINITIONS]
        
        # Compact multi-select for motif classes
        selected_classes = st.multiselect(
            "Choose motif classes to analyze:",
            options=list(CLASS_DEFINITIONS.keys()),
            default=st.session_state.get("selected_classes", []),
            format_func=lambda x: f"{CLASS_DEFINITIONS[x]['name']} ({len(CLASS_DEFINITIONS[x]['subclasses'])} subclasses)",
            help="Select one or more motif classes for analysis"
        )
        st.session_state["selected_classes"] = selected_classes
        
        # Show selection summary
        if selected_classes:
            total_subclasses = sum(len(CLASS_DEFINITIONS[cls]["subclasses"]) for cls in selected_classes)
            st.success(f"✅ {len(selected_classes)} classes selected ({total_subclasses} subclasses)")
        else:
            st.warning("⚠️ No classes selected")
        
        
        # ========== RUN ANALYSIS BUTTON ========== 
        if st.button("🔬 Run Motif Analysis", type="primary", use_container_width=True, key="run_motif_analysis_main"):
            # Validate selection
            if not st.session_state.selected_classes:
                st.error("❌ Please select at least one motif class to analyze.")
                st.session_state.analysis_status = "Error"
            elif not st.session_state.seqs:
                st.error("❌ Please upload or input sequences before running analysis.")
                st.session_state.analysis_status = "Error"
            else:
                st.session_state.analysis_status = "Running"
                
                # Set analysis parameters based on requirements
                nonoverlap = True  # Keep overlaps disabled for specificity
                report_hotspots = True  # Enable hotspot detection 
                calculate_conservation = False  # Disable to reduce computation time
                threshold = 0.0  # Show all detected motifs (even 0 scores)
                
                validation_messages = []

                # Scientific validation check
                if CONFIG_AVAILABLE and st.session_state.get('selected_classes'):
                    for class_id in st.session_state.selected_classes:
                        limits = get_motif_limits(class_id)
                        if limits:
                            validation_messages.append(f"✓ {class_id}: Length limits {limits}")
                
                st.info("⏳ Running analysis...")
                
                try:
                    # Filter which classes to analyze based on selection
                    analysis_classes = st.session_state.selected_classes if st.session_state.selected_classes else None
                    
                    # Run analysis on each sequence
                    all_results = []
                    all_hotspots = []
                    
                    with st.progress(0, text="Analyzing sequences..."):
                        for i, (seq, name) in enumerate(zip(st.session_state.seqs, st.session_state.names)):
                            progress = (i + 1) / len(st.session_state.seqs)
                            
                            # Run the core analysis
                            results = all_motifs_refactored(
                                seq, name,
                                nonoverlap=nonoverlap,
                                report_hotspots=report_hotspots,
                                calculate_conservation=calculate_conservation
                            )
                            
                            # Filter results by selected classes if specified
                            if analysis_classes:
                                filtered_results = []
                                for motif in results:
                                    motif_class = motif.get('Class', motif.get('Type', 'Unknown'))
                                    if motif_class in analysis_classes:
                                        filtered_results.append(motif)
                                results = filtered_results
                            
                            # Ensure all motifs have required fields
                            results = [ensure_subclass(motif) for motif in results]
                            all_results.append(results)
                            
                            # Handle hotspots if present
                            if isinstance(results, tuple) and len(results) == 2:
                                motifs, hotspots = results
                                all_results[-1] = motifs
                                all_hotspots.extend(hotspots)
                            
                            st.progress(progress, text=f"Analyzed {i+1}/{len(st.session_state.seqs)} sequences")
                    
                    # Store results
                    st.session_state.results = all_results
                    st.session_state.hotspots = all_hotspots
                    
                    # Generate summary
                    summary = []
                    for i, results in enumerate(all_results):
                        stats = get_basic_stats(st.session_state.seqs[i])
                        summary.append({
                            'Sequence': st.session_state.names[i],
                            'Length': stats['Length'],
                            'GC Content': f"{stats['GC%']:.1f}%",
                            'Motifs Found': len(results),
                            'Unique Types': len(set(m.get('Type', 'Unknown') for m in results)),
                            'Avg Score': f"{np.mean([m.get('Score', 0) for m in results]):.3f}" if results else "0.000"
                        })
                    
                    st.session_state.summary_df = pd.DataFrame(summary)
                    st.success("✅ Analysis complete! Results are available below and in the 'Analysis Results and Visualization' tab.")
                    st.session_state.analysis_status = "Complete"
                    
                except Exception as e:
                    st.error(f"❌ Analysis failed: {str(e)}")
                    st.session_state.analysis_status = "Error"

        # Show quick summary table if available
        if st.session_state.get('summary_df') is not None:
            st.markdown("#### Analysis Summary")
            st.dataframe(st.session_state.summary_df)
    # End of Upload & Analyze tab
    st.markdown("---")
# ---------- RESULTS ----------
with tab_pages["Results"]:
    st.markdown('<h2>Analysis Results and Visualization</h2>', unsafe_allow_html=True)
    if not st.session_state.results:
        st.info("No analysis results. Please run motif analysis first.")
    else:
        # Enhanced summary display
        st.markdown("### 📊 Analysis Summary")
        st.dataframe(st.session_state.summary_df, use_container_width=True)
        
        # Sequence selection for detailed analysis
        seq_idx = 0
        if len(st.session_state.seqs) > 1:
            seq_idx = st.selectbox("Choose Sequence for Details:", range(len(st.session_state.seqs)), 
                                 format_func=lambda i: st.session_state.names[i])
        
        motifs = st.session_state.results[seq_idx]
        sequence_length = len(st.session_state.seqs[seq_idx])
        sequence_name = st.session_state.names[seq_idx]
        
        if not motifs:
            st.warning("No motifs detected for this sequence.")
        else:
            # Create enhanced motifs DataFrame
            df = pd.DataFrame(motifs)
            
            # Calculate and display enhanced coverage statistics
            stats = get_basic_stats(st.session_state.seqs[seq_idx], motifs)
            
            # Filter motifs for density calculation (exclude hybrid and cluster)
            filtered_motifs = [m for m in motifs if m.get('Class') not in ['Hybrid', 'Non-B_DNA_Cluster']]
            excluded_motifs = [m for m in motifs if m.get('Class') in ['Hybrid', 'Non-B_DNA_Cluster']]
            
            motif_count = len(motifs)
            filtered_count = len(filtered_motifs)
            excluded_count = len(excluded_motifs)
            coverage_pct = stats.get("Motif Coverage %", 0)
            non_b_density = (filtered_count / sequence_length * 1000) if sequence_length > 0 else 0
            
            # Enhanced summary card
            st.markdown(f"""
            <div style='background: linear-gradient(90deg, #667eea 0%, #764ba2 100%); 
                        border-radius: 12px; padding: 15px; margin: 10px 0; color: white;'>
                <h3 style='margin: 0; color: white; text-align: center;'>Comprehensive Non-B DNA Analysis</h3>
                <div style='display: flex; justify-content: space-around; margin-top: 15px; flex-wrap: wrap;'>
                    <div style='text-align: center; min-width: 120px;'>
                        <h2 style='margin: 5px; color: #FFD700;'>{coverage_pct:.2f}%</h2>
                        <p style='margin: 0; font-size: 16px;'>Motif Coverage*</p>
                    </div>
                    <div style='text-align: center; min-width: 120px;'>
                        <h2 style='margin: 5px; color: #FFD700;'>{non_b_density:.2f}</h2>
                        <p style='margin: 0; font-size: 16px;'>Non-B DNA Density*<br>(motifs/kb)</p>
                    </div>
                    <div style='text-align: center; min-width: 120px;'>
                        <h2 style='margin: 5px; color: #FFD700;'>{motif_count}</h2>
                        <p style='margin: 0; font-size: 16px;'>Total Motifs</p>
                    </div>
                    <div style='text-align: center; min-width: 120px;'>
                        <h2 style='margin: 5px; color: #FFD700;'>{sequence_length:,}</h2>
                        <p style='margin: 0; font-size: 16px;'>Sequence Length (bp)</p>
                    </div>
                </div>
               
            </div>
            """, unsafe_allow_html=True)
            
            # Score analysis section (using raw scores as requested)
            if motifs:
                st.markdown("### 📈 Score Analysis")
                
                score_col1, score_col2 = st.columns(2)
                
                with score_col1:
                    # Raw score distribution
                    raw_scores = []
                    
                    for m in motifs:
                        # Get actual/raw score from various possible keys
                        raw_score = m.get('Actual_Score')
                        if raw_score is None:
                            raw_score = m.get('Score', 0)
                        
                        # Convert to float and add if valid
                        try:
                            raw_float = float(raw_score) if raw_score is not None else 0.0
                            raw_scores.append(raw_float)
                        except (ValueError, TypeError):
                            raw_scores.append(0.0)
                    
                    # Raw score distribution analysis
                    if raw_scores and any(s > 0 for s in raw_scores):
                        fig, ax = plt.subplots(figsize=(10, 4))
                        
                        ax.hist(raw_scores, bins=20, alpha=0.7, color='skyblue', label='Raw Scores')
                        ax.set_xlabel('Raw Score')
                        ax.set_ylabel('Frequency')
                        ax.set_title('Raw Score Distribution')
                        ax.legend()
                        
                        plt.tight_layout()
                        st.pyplot(fig)
                        plt.close(fig)
                    else:
                        st.warning("⚠️ All scores are zero. This might indicate an issue with motif scoring.")
                        st.info("This could be due to:\n- Low sequence quality\n- Missing motif features\n- Scoring threshold too high")
                
                with score_col2:
                    # Motif class distribution - show all classes even if 0
                    motif_classes = [m.get('Class', 'Unknown') for m in motifs]
                    
                    # Initialize all classes with 0 counts  
                    all_class_counts = {cls: 0 for cls in MOTIF_ORDER}
                    
                    # Update with actual counts
                    actual_counts = Counter(motif_classes)
                    all_class_counts.update(actual_counts)
                    
                    # Remove 'Unknown' if it exists and has 0 count
                    if 'Unknown' in all_class_counts and all_class_counts['Unknown'] == 0:
                        del all_class_counts['Unknown']
                    
                    fig, ax = plt.subplots(figsize=(10, 8))
                    
                    # Create labels showing count and percentage
                    labels = []
                    values = list(all_class_counts.values())
                    total = sum(values) if sum(values) > 0 else 1  # Avoid division by zero
                    
                    for cls, count in all_class_counts.items():
                        percentage = (count / total) * 100 if total > 0 else 0
                        if count > 0:
                            labels.append(f'{cls}\n({count}, {percentage:.1f}%)')
                        else:
                            labels.append(f'{cls}\n(0, 0.0%)')
                    
                    # Use colors from MOTIF_COLORS where available
                    colors = []
                    for cls in all_class_counts.keys():
                        if cls in MOTIF_COLORS:
                            colors.append(MOTIF_COLORS[cls])
                        else:
                            # Default colors for any missing classes
                            default_colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99', 
                                            '#ff99cc', '#99ccff', '#ffccbb', '#ccffcc']
                            colors.append(default_colors[len(colors) % len(default_colors)])
                    
                    # Create bar chart instead of pie chart for better readability when many classes have 0
                    bars = ax.bar(range(len(all_class_counts)), values, color=colors)
                    ax.set_xlabel('Motif Classes')
                    ax.set_ylabel('Count')
                    ax.set_title('Motif Class Distribution (All Classes)')
                    ax.set_xticks(range(len(all_class_counts)))
                    ax.set_xticklabels(list(all_class_counts.keys()), rotation=45, ha='right')
                    
                    # Add count labels on bars
                    for bar, count in zip(bars, values):
                        if count > 0:
                            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
                                   str(count), ha='center', va='bottom')
                    
                    plt.tight_layout()
                    st.pyplot(fig)
                    plt.close(fig)
            
            # Enhanced motif table with new columns
            st.markdown(f"### 📋 Detailed Motif Table for **{sequence_name}**")
            
            # Column selection for display
            available_columns = df.columns.tolist()
            display_columns = st.multiselect(
                "Select columns to display:",
                available_columns,
                default=[col for col in ['Class', 'Subclass', 'Start', 'End', 'Length', 'Actual Score', 'Score', 'GC Content'] if col in available_columns]
            )
            
            if display_columns:
                # Filter out sequence and sequence name columns for cleaner display
                filtered_df = df[display_columns].copy()
                # Replace underscores with spaces in column names for display
                filtered_df.columns = [col.replace('_', ' ') for col in filtered_df.columns]
                st.dataframe(filtered_df, use_container_width=True, height=360)
            else:
                # Replace underscores with spaces in column names for display
                display_df = df.copy()
                display_df.columns = [col.replace('_', ' ') for col in display_df.columns]
                st.dataframe(display_df, use_container_width=True, height=360)
            
            # AUTOMATIC COMPREHENSIVE VISUALIZATION GENERATION
            st.markdown('<h3>📊 Comprehensive Analysis - Information-Based Visualizations</h3>', unsafe_allow_html=True)
            st.info("🎯 Generating all visualizations automatically (organized by information type, not plot type)")
            
            with st.spinner("Creating comprehensive information-based visualization suite..."):
                try:
                    # Generate comprehensive visualizations automatically
                    static_plots, interactive_plots, detailed_stats = create_comprehensive_information_based_visualizations(
                        df, sequence_length, sequence_name)
                    
                    # Display information-type organized visualizations
                    visualization_categories = [
                        ("📈 Coverage & Density Analysis", ["coverage_analysis", "detailed_coverage_map"]),
                        ("📊 Distribution Analysis", ["distribution_analysis"]),
                        ("🧬 Sequence Analysis", ["sequence_analysis"]),
                        ("⚖️ Comparative Analysis", ["comparative_analysis"]),
                        ("🔬 Advanced Analysis", ["advanced_analysis"])
                    ]
                    
                    for category_name, plot_keys in visualization_categories:
                        st.markdown(f'<h4>{category_name}</h4>', unsafe_allow_html=True)
                        
                        # Display plots in this category
                        for plot_key in plot_keys:
                            if plot_key in static_plots:
                                st.pyplot(static_plots[plot_key])
                                plt.close(static_plots[plot_key])  # Free memory
                    
                    # Interactive Visualizations Section
                    if interactive_plots:
                        st.markdown('<h4>🎯 Interactive Visualizations</h4>', unsafe_allow_html=True)
                        
                        for plot_name, plot_fig in interactive_plots.items():
                            st.plotly_chart(plot_fig, use_container_width=True)
                    
                    st.success("✅ All visualizations generated successfully! All plots are organized by information type.")
                    
                except Exception as e:
                    st.error(f"Error generating comprehensive visualizations: {e}")
                    st.info("Falling back to basic visualization options...")
                    
                    # Fallback to basic visualization selector
                    viz_options = [
                        "Basic Charts", "Interactive Plots", "Statistical Analysis", 
                        "Genomic Mapping", "Advanced Analysis"
                    ]
                    selected_viz = st.radio("Choose Visualization Category:", viz_options, horizontal=True)
            
                    # Simple fallback visualization for error cases
                    st.markdown('<h4>📊 Fallback Visualization</h4>', unsafe_allow_html=True)
                    st.info("Using basic fallback visualizations since comprehensive analysis failed.")
                    
                    # Basic motif count chart
                    st.markdown('<span style="font-family:Montserrat,Arial;font-size:1.11rem;"><b>Motif Class Distribution</b></span>', unsafe_allow_html=True)
                    fig, ax = plt.subplots(figsize=(10,6))
                    class_counts = df['Class'].value_counts()
                    ax.barh(class_counts.index, class_counts.values)
                    ax.set_xlabel("Motif Count")
                    ax.set_title("Basic Motif Class Distribution")
                    st.pyplot(fig)
                    plt.close(fig)
                    
                    # Basic motif map
                    st.markdown('<span style="font-family:Montserrat,Arial;font-size:1.11rem;"><b>Motif Position Map</b></span>', unsafe_allow_html=True)
                    fig, ax = plt.subplots(figsize=(12,4))
                    for i, (_, row) in enumerate(df.iterrows()):
                        ax.plot([row['Start'], row['End']], [i, i], lw=3, alpha=0.7)
                    ax.set_xlabel("Sequence Position (bp)")
                    ax.set_ylabel("Motif Index")
                    ax.set_title("Basic Motif Position Map")
                    st.pyplot(fig)
                    plt.close(fig)

# ---------- DOWNLOAD ----------
with tab_pages["Download"]:
    st.header("Export Data")
    if not st.session_state.results:
        st.info("No results available to download.")
    else:
        st.markdown("### 📊 Export Options")
        
        # Export configuration
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**📊 Export Configuration**")
            
        with col2:
            include_sequences = st.checkbox("Include Full Sequences", value=True,
                                           help="Include full motif sequences in export")
        
        # Prepare data for export
        df_all = []
        for i, motifs in enumerate(st.session_state.results):
            for m in motifs:
                export_row = m.copy()
                # Only add Sequence Name if it's not already present
                if 'Sequence Name' not in export_row:
                    export_row['Sequence Name'] = st.session_state.names[i]
                
                # Handle class mapping for compatibility
                if export_row['Class'] == "Z-DNA" and export_row.get("Subclass", "") == "eGZ (Extruded-G)":
                    export_row['Class'] = "eGZ (Extruded-G)"
                
                # Filter sequence data if not requested
                if not include_sequences:
                    export_row.pop('Sequence', None)
                
                # Keep both score types (normalized and actual)
                # No filtering of score columns
                
                df_all.append(export_row)
        
        df_all = pd.DataFrame(df_all)
        
        # Remove any duplicate columns that might have been created
        df_all = df_all.loc[:, ~df_all.columns.duplicated()]
        
        # Replace underscores with spaces in column names for better readability
        df_all.columns = [col.replace('_', ' ') for col in df_all.columns]
        
        # Remove any duplicate columns that might have been created after underscore replacement
        df_all = df_all.loc[:, ~df_all.columns.duplicated()]
        
        # Display preview
        st.markdown("### 👀 Export Preview")
        st.dataframe(df_all.head(10), use_container_width=True)
        st.caption(f"Showing first 10 of {len(df_all)} total records")
        
        # Export buttons - Always provide both CSV and Excel
        st.markdown("### 💾 Download Files")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            csv_data = df_all.to_csv(index=False).encode("utf-8")
            st.download_button(
                "📄 Download CSV", 
                data=csv_data, 
                file_name="nbdfinder_results_both_scores.csv", 
                mime="text/csv",
                use_container_width=True,
                help="CSV format with both normalized and actual scores"
            )
        
        with col2:
            excel_data = io.BytesIO()
            with pd.ExcelWriter(excel_data, engine='xlsxwriter') as writer:
                df_all.to_excel(writer, index=False, sheet_name="Motifs")
            
            excel_data.seek(0)
            st.download_button(
                "📊 Download Excel", 
                data=excel_data, 
                file_name="nbdfinder_results_both_scores.xlsx", 
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                use_container_width=True,
                help="Excel format with both normalized and actual scores"
            )
        
        with col3:
            # Enhanced export options with new formats
            if CONFIG_AVAILABLE:
                config_summary = {
                    'Analysis Parameters': {
                        'Score Type': 'Both (Normalized and Actual)',
                        'Normalized Scoring': True,  # Always use normalized scoring
                        'Score Threshold': 0.0,  # Include all detected motifs
                        'Overlaps Removed': 'Within class only',  # Overlaps prevented within class, allowed between classes
                        'Hotspots Detected': True  # Always detect hotspots
                    },
                    'Motif Length Limits': MOTIF_LENGTH_LIMITS,
                    'Scoring Methods': SCORING_METHODS
                }
                
                import json
                config_json = json.dumps(config_summary, indent=2)
                st.download_button(
                    "⚙️ Download Config", 
                    data=config_json, 
                    file_name="analysis_configuration.json", 
                    mime="application/json",
                    use_container_width=True
                )
        
        # Prepare data for genome browser exports
        all_seq_data = []
        for i, motifs in enumerate(st.session_state.results):
            all_seq_data.append({
                'sequence_name': st.session_state.names[i],
                'motifs': motifs
            })
        
        # NEW: Genome Browser Export Section
        st.markdown("### 🧬 Genome Browser Formats")
        st.caption("Export data for UCSC Genome Browser, IGV, and other genomic tools")
        
        col4, col5, col6 = st.columns(3)
        
        with col4:
            # BED format export
            try:
                from nbdio.writers import export_to_bed
                bed_data = export_to_bed(
                    [motif for seq_data in all_seq_data for motif in seq_data['motifs']], 
                    sequence_name="NBDFinder_Analysis",
                    score_type="normalized",
                    include_subclass=True
                )
                st.download_button(
                    "📄 Download BED",
                    data=bed_data,
                    file_name="nbdfinder_motifs_normalized.bed",
                    mime="text/plain",
                    use_container_width=True,
                    help="BED format for UCSC Genome Browser and IGV (normalized scores)"
                )
            except Exception as e:
                st.error(f"BED export error: {e}")
        
        with col5:
            # GFF3 format export  
            try:
                from nbdio.writers import export_to_gff3
                gff3_data = export_to_gff3(
                    [motif for seq_data in all_seq_data for motif in seq_data['motifs']], 
                    sequence_name="NBDFinder_Analysis"
                )
                st.download_button(
                    "📋 Download GFF3",
                    data=gff3_data,
                    file_name="nbdfinder_motifs.gff3",
                    mime="text/plain",
                    use_container_width=True,
                    help="GFF3 format for detailed genomic annotations"
                )
            except Exception as e:
                st.error(f"GFF3 export error: {e}")
        
        with col6:
            # Density bedGraph export
            try:
                from nbdio.writers import create_density_bedgraph
                if all_seq_data:
                    # Filter out hybrid and cluster motifs for density export
                    all_motifs = [motif for seq_data in all_seq_data for motif in seq_data['motifs']]
                    filtered_motifs = [m for m in all_motifs if m.get('Class') not in ['Hybrid', 'Non-B_DNA_Cluster']]
                    
                    # Use the first sequence for length estimation
                    seq_length = max([motif.get('End', 0) for motif in all_motifs], default=1000)
                    density_data = create_density_bedgraph(
                        filtered_motifs,
                        seq_length,
                        sequence_name="NBDFinder_Analysis"
                    )
                    
                    excluded_count = len(all_motifs) - len(filtered_motifs)
                    help_text = "BedGraph format for motif density visualization"
                    if excluded_count > 0:
                        help_text += f" (excludes {excluded_count} hybrid/cluster motifs)"
                    
                    st.download_button(
                        "📊 Download Density",
                        data=density_data,
                        file_name="nbdfinder_density.bedgraph",
                        mime="text/plain", 
                        use_container_width=True,
                        help=help_text
                    )
            except Exception as e:
                st.error(f"Density export error: {e}")
                
        # Cache statistics display
        st.markdown("### 📊 Performance & Caching")
        try:
            from enhanced_cache import get_cache_stats
            cache_stats = get_cache_stats()
            
            col7, col8, col9, col10 = st.columns(4)
            with col7:
                st.metric("Cache Hit Rate", f"{cache_stats['totals']['hit_rate_percent']:.1f}%")
            with col8:
                st.metric("Memory Used", f"{cache_stats['totals']['memory_used_mb']:.1f} MB")
            with col9:
                st.metric("Total Requests", cache_stats['totals']['total_requests'])
            with col10:
                if st.button("🗑️ Clear Cache", use_container_width=True):
                    from enhanced_cache import clear_all_caches
                    clear_all_caches()
                    st.success("All caches cleared!")
                    st.experimental_rerun()
                    
        except ImportError:
            st.info("Enhanced caching not available")

# ---------- DOCUMENTATION ----------
with tab_pages["Documentation"]:
    st.header("Scientific Documentation & References")
    
    # Add new visualization documentation
    st.markdown("""
    <div style='background:#f4faff; border-radius:12px; padding:18px; font-size:1.08rem; font-family:Montserrat,Arial;'>
    <b>🎨 Enhanced Visualization Suite (New Features)</b><br><br>
    
    The NBDFinder tool now includes a comprehensive visualization suite with 21+ chart types organized into 5 categories:
    
    <ul>
        <li><b>Basic Charts</b>: Motif counts, pie charts, stacked distributions, and motif tracks</li>
        <li><b>Interactive Plots</b>: Plotly-powered sunburst, interactive browsers, and track plots</li>
        <li><b>Statistical Analysis</b>: Score distributions, CDF plots, t-SNE clustering, and Manhattan plots</li>
        <li><b>Genomic Mapping</b>: Position analysis, density heatmaps, sequence coverage, and GC content scatter</li>
        <li><b>Advanced Analysis</b>: Class-subclass heatmaps, network graphs, Venn diagrams, and cluster density</li>
    </ul>
    
    <b>🔧 Recent Updates & Integration</b><br>
    These visualization features were developed in recent PRs but were previously not integrated into the Streamlit interface. 
    They are now fully accessible through the "Analysis Results and Visualization" tab with an intuitive category-based selector.
    </div>
    <br>
    """, unsafe_allow_html=True)
    
    # NEW: REST API Documentation
    st.markdown("""
    <div style='background:#fff4e6; border-radius:12px; padding:18px; font-size:1.08rem; font-family:Montserrat,Arial;'>
    <b>🚀 REST API Access (New Feature)</b><br><br>
    
    NBDFinder now provides a RESTful API for programmatic access and integration with other tools:
    
    <ul>
        <li><b>Start API Server</b>: <code>python api.py</code> (runs on http://localhost:8000)</li>
        <li><b>Health Check</b>: <code>GET /api/v1/health</code></li>
        <li><b>Analyze Sequence</b>: <code>POST /api/v1/analyze</code></li>
        <li><b>Get Statistics</b>: <code>GET /api/v1/stats</code></li>
        <li><b>Motif Information</b>: <code>GET /api/v1/motif-info</code></li>
    </ul>
    
    <b>Example Usage:</b><br>
    <code>
    curl -X POST "http://localhost:8000/api/v1/analyze" \\<br>
    &nbsp;&nbsp;-H "Content-Type: application/json" \\<br>
    &nbsp;&nbsp;-d '{"sequence": "GGGTTAGGGTTAGGGTTAGGG", "sequence_name": "test"}'
    </code><br><br>
    
    <b>Features:</b> JSON responses, comprehensive motif data, caching, CORS support, and automatic API documentation at <code>/docs</code>
    </div>
    <br>
    """, unsafe_allow_html=True)
    
    st.markdown("""
    <div style='background:#f4faff; border-radius:12px; padding:18px; font-size:1.08rem; font-family:Montserrat,Arial;'>
    <b>Motif Classes Detected:</b><br><br>
    <ul>
        <li><b>Curved DNA</b>: Identifies phased poly(A) or poly(T) tracts using regex and spacing rules, reflecting intrinsic curvature. Scoring is based on tract length/grouping.</li>
        <li><b>Z-DNA</b>: Detects alternating purine-pyrimidine patterns, GC-rich segments. Uses windowed scoring; regex finds dinucleotide repeats.</li>
        <li><b>eGZ-motif (Extruded-G Z-DNA)</b>: Searches for long (CGG)<sub>n</sub> runs via regex. Scored by repeat count.</li>
        <li><b>Slipped DNA</b>: Recognizes direct/tandem repeats by repeat-unit matching and regex. Scoring by length and unit copies.</li>
        <li><b>R-Loop</b>: Finds G-rich regions for stable RNA-DNA hybrids; RLFS model and regex. Thermodynamic scoring for hybrid stability.</li>
        <li><b>Cruciform</b>: Finds palindromic inverted repeats with spacers, regex and reverse complement. Scoring by arm length and A/T content.</li>
        <li><b>Triplex DNA / Mirror Repeat</b>: Detects purine/pyrimidine mirror repeats/triplex motifs. Regex identifies units; scoring by composition/purity.</li>
        <li><b>Sticky DNA</b>: Searches extended GAA/TTC repeats. Scoring by repeat count.</li>
        <li><b>G-Triplex</b>: Finds three consecutive guanine runs by regex and loop length. Scoring by G-run sum and loop penalty.</li>
        <li><b>G4 (G-Quadruplex) and Variants</b>: Detects canonical/variant G4 motifs by G-run/loop regex. G4Hunter scoring for content/structure.</li>
        <li><b>i-Motif</b>: C-rich sequences for i-motif under acid. Regex for C runs/loops; scoring by run count and content.</li>
        <li><b>AC-Motif</b>: Alternating A-rich/C-rich consensus regions by regex. Scoring by pattern presence.</li>
        <li><b>A-philic DNA</b>: Uses tetranucleotide log2 odds scoring to identify A-tract-favoring sequences with high protein-binding affinity. Classified as high-confidence or moderate A-philic based on score thresholds and strong tetranucleotide counts.</li>
        <li><b>Hybrid Motif</b>: Regions where motif classes overlap; found by interval intersection, scored on diversity/size.</li>
        <li><b>Non-B DNA Clusters</b>: Hotspots with multiple motifs in a window; sliding algorithm, scored by motif count/diversity.</li>
    </ul>
    <b>References:</b>
    <ul>
        <li>Bedrat et al., 2016 Nucleic Acids Research</li>
        <li>Ho et al., 2010 Nature Chemical Biology</li>
        <li>Kim et al., 2018 Nucleic Acids Research</li>
        <li>Zeraati et al., 2018 Nature Chemistry</li>
        <li>Bacolla et al., 2006 Nucleic Acids Research</li>
        <li>Mirkin & Frank-Kamenetskii, 1994 Annual Review of Biophysics</li>
        <li>Vinogradov, 2003 Bioinformatics (A-philic DNA tetranucleotide analysis)</li>
        <li>Bolshoy et al., 1991 PNAS (A-tract structural properties)</li>
        <li>Rohs et al., 2009 Nature (Protein-DNA interactions, A-philic binding)</li>
        <li>New et al., 2020 Journal of DNA Structure</li>
    </ul>
    </div>
    """, unsafe_allow_html=True)
    
    # Add configuration information if available
    if CONFIG_AVAILABLE:
        st.markdown("""
        <div style='background:#f1f5f9; border-radius:12px; padding:18px; font-size:1.08rem; font-family:Montserrat,Arial; margin-top:20px;'>
        <b>📋 Scoring Configuration Details</b>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown("### Motif Length Constraints")
        
        config_df = pd.DataFrame([
            {
                "Motif Class": motif_class,
                "Min Length (bp)": limits["S_min"],
                "Max Length (bp)": limits["S_max"],
                "Biological Basis": f"Based on {SCORING_METHODS.get(motif_class, {}).get('reference', 'literature survey')}"
            }
            for motif_class, limits in MOTIF_LENGTH_LIMITS.items()
        ])
        
        st.dataframe(config_df, use_container_width=True)
        
        st.markdown("### Scoring Methods")
        scoring_df = pd.DataFrame([
            {
                "Motif Class": motif_class,
                "Method": method_info.get("method", ""),
                "Description": method_info.get("description", ""),
                "Reference": method_info.get("reference", "")
            }
            for motif_class, method_info in SCORING_METHODS.items()
        ])
        
        st.dataframe(scoring_df, use_container_width=True)

st.markdown("""
---
<div style='font-size: 1.05rem; color: #1e293b; margin-top: 30px; text-align: left; font-family:Montserrat,Arial;'>
<b>Developed by</b><br>
Dr. Venkata Rajesh Yella<br>
<a href='mailto:yvrajesh_bt@kluniversity.in'>yvrajesh_bt@kluniversity.in</a> |
<a href='https://github.com/VRYella' target='_blank'>GitHub: VRYella</a>
</div>
""", unsafe_allow_html=True)
