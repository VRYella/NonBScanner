# --- Imports and Config ---
import streamlit as st, pandas as pd, matplotlib.pyplot as plt, io, numpy as np
from collections import Counter; from Bio import Entrez, SeqIO
from utils import parse_fasta, gc_content, reverse_complement, wrap
from motifs import get_basic_stats
from motifs.base_motif import select_best_nonoverlapping_motifs
from all_motifs_refactored import all_motifs_refactored
from motifs.enhanced_visualization import create_comprehensive_information_based_visualizations

# --- Optional Config Import ---
try:
    from classification_config import MOTIF_LENGTH_LIMITS, SCORING_METHODS, get_motif_limits; CONFIG_AVAILABLE = True
except ImportError:
    CONFIG_AVAILABLE = False; MOTIF_LENGTH_LIMITS = {}; SCORING_METHODS = {}

# --- Streamlit Page Setup ---
st.set_page_config(page_title="Non-B DNA Motif Finder", layout="wide", page_icon="🧬",
    menu_items={'About': "Non-B DNA Motif Finder | Developed by Dr. Venkata Rajesh Yella"})

# --- Motif subclass patch ---
def ensure_subclass(m):
    if isinstance(m, dict) and ('Subclass' not in m or m['Subclass'] is None):
        m['Subclass'] = m.get('Subtype', 'Other'); return m
    return m if isinstance(m, dict) else {'Subclass': 'Other', 'Motif': m}

# --- CSS ---
st.markdown("""<style>
body,[data-testid="stAppViewContainer"],.main{background:#f7fafd!important;font-family:'Montserrat',Arial,sans-serif!important;}
.stTabs [data-baseweb="tab-list"]{width:100vw!important;justify-content:stretch!important;border-bottom:2px solid #1565c0;background:linear-gradient(90deg,#eaf3fa 0%,#f7fafd 100%)!important;box-shadow:0 2px 8px #dae5f2;}
.stTabs [data-baseweb="tab"]{font-size:1.45rem!important;font-weight:700!important;flex:1 1 0%;padding:15px 0!important;text-align:center;color:#1565c0!important;background:#eaf3fa!important;border-right:1px solid #eee!important;}
.stTabs [aria-selected="true"]{color:#002147!important;border-bottom:5px solid #1565c0!important;background:#f7fafd!important;box-shadow:0 4px 8px #e0e5ea;}
.stTabs [data-baseweb="tab"]:last-child{border-right:none!important;}
h1,h2,h3,h4{font-family:'Montserrat',Arial,sans-serif!important;color:#1565c0!important;font-weight:800!important;}
h1{font-size:2.05rem!important;}h2{font-size:1.55rem!important;}h3{font-size:1.19rem!important;}h4{font-size:1.09rem!important;}
.stMarkdown,.markdown-text-container,.stText,p,span,label,input,.stTextInput>div>div>input,.stSelectbox>div>div>div,.stMultiSelect>div>div>div,.stRadio>div>div>label>div{font-size:1.08rem!important;font-family:'Montserrat',Arial,sans-serif!important;line-height:1.6!important;}
.stButton>button{font-size:1.08rem!important;font-family:'Montserrat',Arial,sans-serif!important;padding:0.45em 1.2em!important;background:linear-gradient(90deg,#1565c0 0%,#2e8bda 100%)!important;color:#fff!important;border-radius:7px!important;border:none!important;font-weight:600!important;}
.stButton>button:hover{background:linear-gradient(90deg,#2e8bda 0%,#1565c0 100%)!important;}
.stDataFrame,.stTable{font-size:1.05rem!important;font-family:'Montserrat',Arial,sans-serif!important;}
.stTabs [data-baseweb="tab-panel"]{padding-top:2rem!important;padding-left:1rem!important;padding-right:1rem!important;}
.analysis-summary-card{margin:1.5rem 0!important;padding:1.5rem!important;}
.stSelectbox > div > div > div,.stMultiSelect > div > div > div{min-height:3.0rem!important;padding:0.8rem!important;line-height:1.5!important;}
.stTextInput > div > div > input,.stTextArea textarea{padding:0.75rem 1rem!important;min-height:2.8rem!important;line-height:1.4!important;}
.stRadio > div > div > label{margin-bottom:0.5rem!important;padding:0.3rem 0.5rem!important;}
.stNumberInput > div > div > input{padding:0.75rem 1rem!important;min-height:2.8rem!important;}
</style>""", unsafe_allow_html=True)

# --- App Constants & Motif Order/Colors ---
if CONFIG_AVAILABLE:
    expanded_order = []
    for m in list(MOTIF_LENGTH_LIMITS.keys()):
        if m == "Slipped_DNA_DR": expanded_order += ["Slipped DNA (Direct Repeat)", "Slipped DNA (STR)"]
        elif m == "Slipped_DNA_STR": continue
        elif m == "eGZ": expanded_order.append("eGZ (Extruded-G)")
        elif m == "G4": expanded_order += ["G4", "Relaxed G4", "Bulged G4", "Bipartite G4", "Multimeric G4"]
        elif m == "G-Triplex": expanded_order.append("G-Triplex")
        elif m == "AC-motif": expanded_order.append("AC-Motif")
        elif m == "A-philic_DNA": expanded_order.append("A-philic DNA")
        else:
            disp = m.replace("_"," ").replace("-","-")
            if m=="Curved_DNA": disp="Curved DNA"
            elif m=="Z-DNA": disp="Z-DNA"
            elif m=="R-Loop": disp="R-Loop"
            elif m=="Triplex": disp="Triplex DNA"
            elif m=="Sticky_DNA": disp="Sticky DNA"
            elif m=="i-Motif": disp="i-Motif"
            expanded_order.append(disp)
    MOTIF_ORDER = expanded_order+["Hybrid","Non-B DNA Clusters"]
else:
    MOTIF_ORDER = ["Sticky DNA","Curved DNA","Z-DNA","eGZ (Extruded-G)","Slipped DNA","R-Loop","Cruciform","Triplex DNA","G-Triplex","G4","Relaxed G4","Bulged G4","Bipartite G4","Multimeric G4","i-Motif","AC-Motif","A-philic DNA","Hybrid","Non-B DNA Clusters"]

MOTIF_COLORS = {"Curved DNA":"#FF9AA2","Z-DNA":"#FFB7B2","eGZ (Extruded-G)":"#6A4C93","Slipped DNA":"#FFDAC1","Slipped DNA (Direct Repeat)":"#FFDAC1","Slipped DNA (STR)":"#FFE4B3","R-Loop":"#FFD3B6","Cruciform":"#E2F0CB","Triplex DNA":"#B5EAD7","Sticky DNA":"#DCB8CB","G-Triplex":"#C7CEEA","G4":"#A2D7D8","Relaxed G4":"#A2D7B8","Bulged G4":"#A2A7D8","Bipartite G4":"#A2D788","Multimeric G4":"#A2A7B8","i-Motif":"#B0C4DE","A-philic DNA":"#E6B8F7","Hybrid":"#C1A192","Non-B DNA Clusters":"#A2C8CC","AC-Motif":"#F5B041"}
PAGES = {"Home":"Overview","Upload & Analyze":"Sequence Upload and Motif Analysis","Results":"Analysis Results and Visualization","Download":"Export Data","Documentation":"Scientific Documentation & References"}
Entrez.email = "raazbiochem@gmail.com"; Entrez.api_key = None
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

# --- Session state setup ---
for k, v in {'seqs':[],'names':[],'results':[],'summary_df':pd.DataFrame(),'hotspots':[],'analysis_status':"Ready"}.items():
    if k not in st.session_state: st.session_state[k]=v

# --- Tabs ---
tabs=st.tabs(list(PAGES.keys())); tab_pages=dict(zip(PAGES.keys(),tabs))

# --- Home Tab ---
with tab_pages["Home"]:
    st.markdown("<h1>Non-B DNA Motif Finder</h1>",unsafe_allow_html=True)
    left,right=st.columns([1,1])
    with left:
        try: st.image("nbdcircle.JPG")
        except: st.markdown("""<div style='background:linear-gradient(135deg,#667eea 0%,#764ba2 100%);border-radius:15px;padding:40px;text-align:center;color:white;'><h2 style='margin:0;color:white;'>🧬</h2><h3 style='margin:10px 0 0 0;color:white;'>NBD Finder</h3><p style='margin:5px 0 0 0;color:#E8E8E8;'>Non-B DNA Detection</p></div>""",unsafe_allow_html=True)
    with right:
        st.markdown("""<div style='font-family:Montserrat,Arial;font-size:1.14rem;color:#222;line-height:1.7;padding:18px;background:#f8f9fa;border-radius:14px;box-shadow:0 2px 8px #eee;'><b>Non-canonical DNA structures</b> play key roles... <b>10 major classes with 22+ subclasses</b> ...</div>""",unsafe_allow_html=True)

# --- Upload & Analyze Tab ---
with tab_pages["Upload & Analyze"]:
    st.markdown("<h2>Sequence Upload and Motif Analysis</h2>",unsafe_allow_html=True)
    st.markdown('<span style="font-family:Montserrat,Arial;font-size:1.12rem;">Supports multi-FASTA and single FASTA. Paste, upload, select example, or fetch from NCBI.</span>',unsafe_allow_html=True)
    st.caption("Supported formats: .fa, .fasta, .txt | Limit: 200MB/file."); st.markdown("---")
    col_left,col_right=st.columns([1,1],gap="large")
    with col_left:
        st.markdown("### 📁 Input Method")
        input_method=st.radio("Choose your input method:",["📂 Upload FASTA File","✏️ Paste Sequence","🧪 Example Data","🌐 NCBI Fetch"],horizontal=True)
        seqs,names=[],[]
        if input_method=="📂 Upload FASTA File":
            fasta_file=st.file_uploader("Drag and drop FASTA/multi-FASTA file here",type=["fa","fasta","txt"])
            if fasta_file:
                content=fasta_file.read().decode("utf-8"); seqs,names=[],[]; cur_seq,cur_name="",""
                for line in content.splitlines():
                    if line.startswith(">"):
                        if cur_seq: seqs.append(parse_fasta(cur_seq)); names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                        cur_name=line.strip().lstrip(">"); cur_seq=""
                    else: cur_seq+=line.strip()
                if cur_seq: seqs.append(parse_fasta(cur_seq)); names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                if seqs:
                    st.success(f"✅ Loaded {len(seqs)} sequences.")
                    for i,seq in enumerate(seqs[:3]):
                        stats=get_basic_stats(seq)
                        st.markdown(f"**{names[i]}**: <span style='color:#576574'>{len(seq):,} bp</span>",unsafe_allow_html=True)
                        st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                    if len(seqs)>3: st.caption(f"...and {len(seqs)-3} more.")
                else: st.warning("No sequences found.")
        elif input_method=="✏️ Paste Sequence":
            seq_input=st.text_area("Paste single or multi-FASTA here:",height=150,placeholder="Paste your DNA sequence(s) here...")
            if seq_input:
                seqs,names=[],[]; cur_seq,cur_name="",""
                for line in seq_input.splitlines():
                    if line.startswith(">"):
                        if cur_seq: seqs.append(parse_fasta(cur_seq)); names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                        cur_name=line.strip().lstrip(">"); cur_seq=""
                    else: cur_seq+=line.strip()
                if cur_seq: seqs.append(parse_fasta(cur_seq)); names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                if seqs:
                    st.success(f"✅ Pasted {len(seqs)} sequences.")
                    for i,seq in enumerate(seqs[:3]):
                        stats=get_basic_stats(seq)
                        st.markdown(f"**{names[i]}**: <span style='color:#576574'>{len(seq):,} bp</span>",unsafe_allow_html=True)
                        st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                    if len(seqs)>3: st.caption(f"...and {len(seqs)-3} more.")
                else: st.warning("No sequences found.")
        elif input_method=="🧪 Example Data":
            ex_type=st.radio("Example Type:",["Single Example","Multi-FASTA Example"],horizontal=True)
            if ex_type=="Single Example" and st.button("🔬 Load Single Example"):
                seqs=[parse_fasta(EXAMPLE_FASTA)]; names=["Example Sequence"]; st.success("✅ Single example sequence loaded.")
                stats=get_basic_stats(seqs[0]); st.code(EXAMPLE_FASTA,language="fasta"); st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
            elif ex_type=="Multi-FASTA Example" and st.button("🔬 Load Multi-FASTA Example"):
                seqs,names=[],[]; cur_seq,cur_name="",""
                for line in EXAMPLE_MULTI_FASTA.splitlines():
                    if line.startswith(">"):
                        if cur_seq: seqs.append(parse_fasta(cur_seq)); names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                        cur_name=line.strip().lstrip(">"); cur_seq=""
                    else: cur_seq+=line.strip()
                if cur_seq: seqs.append(parse_fasta(cur_seq)); names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                st.success(f"✅ Multi-FASTA example loaded with {len(seqs)} sequences.")
                for i,seq in enumerate(seqs[:3]):
                    stats=get_basic_stats(seq)
                    st.markdown(f"**{names[i]}**: <span style='color:#576574'>{len(seq):,} bp</span>",unsafe_allow_html=True)
                    st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                st.code(EXAMPLE_MULTI_FASTA,language="fasta")
        elif input_method=="🌐 NCBI Fetch":
            db=st.radio("NCBI Database",["nucleotide","gene"],horizontal=True)
            st.radio("Query Type",["Accession","Gene Name","Custom Query"],horizontal=True)
            with st.expander("Motif Example Queries"):
                for motif,example in {"G-quadruplex":"NR_003287.2","Z-DNA":"NM_001126112.2","R-loop":"NR_024540.1","eGZ-motif":"CGG repeat region","AC-motif":"A-rich/C-rich consensus region"}.items():
                    st.write(f"**{motif}**: `{example}`")
            query=st.text_input("Enter query (accession, gene, etc.):"); rettype=st.radio("Return Format",["fasta","gb"],horizontal=True)
            retmax=st.number_input("Max Records",min_value=1,max_value=20,value=3)
            if st.button("Fetch from NCBI") and query:
                with st.spinner("Contacting NCBI..."):
                    handle=Entrez.efetch(db=db,id=query,rettype=rettype,retmode="text")
                    records=list(SeqIO.parse(handle,rettype)); handle.close()
                    seqs=[str(rec.seq).upper().replace("U","T") for rec in records]; names=[rec.id for rec in records]
                if seqs:
                    st.success(f"Fetched {len(seqs)} sequences.")
                    for i,seq in enumerate(seqs[:3]):
                        stats=get_basic_stats(seq)
                        st.markdown(f"<b>{names[i]}</b>: <span style='color:#576574'>{len(seq):,} bp</span>",unsafe_allow_html=True)
                        st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                else: st.warning("Enter a query before fetching.")
        if seqs: st.session_state.seqs,st.session_state.names,st.session_state.results=seqs,names,[]
        if st.session_state.get('seqs'):
            st.markdown("### 📊 Sequence Preview")
            for i,seq in enumerate(st.session_state.seqs[:2]):
                stats=get_basic_stats(seq)
                st.markdown(f"**{st.session_state.names[i]}** ({len(seq):,} bp) | GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}",unsafe_allow_html=True)
                st.code(wrap(seq[:400]),language="fasta")
            if len(st.session_state.seqs)>2: st.caption(f"...and {len(st.session_state.seqs)-2} more.")
    with col_right:
        st.markdown("### 🚀 Analysis & Run")
        st.markdown(f"**Configuration Summary:**\n- Motifs: {len(MOTIF_ORDER)} classes selected\n- Scoring: Normalized (0-1)\n- Threshold: 0.0 (all motifs)\n- Overlaps: Not within class\n- Hotspots: Detected")
        nonoverlap,report_hotspots,calc_conservation,threshold=True,True,False,0.0
        if st.button("🔬 Run Motif Analysis",type="primary",use_container_width=True):
            st.session_state.analysis_status="Running"; val_msgs=[]
            for i,seq in enumerate(st.session_state.seqs):
                seq_name=st.session_state.names[i] if i<len(st.session_state.names) else f"Sequence_{i+1}"
                valid_chars=set('ATCGN'); seq_chars=set(seq.upper())
                if not seq_chars.issubset(valid_chars): val_msgs.append(f"⚠️ {seq_name}: Contains non-DNA characters: {seq_chars-valid_chars}")
                if len(seq)<10: val_msgs.append(f"⚠️ {seq_name}: Sequence too short (<10 bp)")
                elif len(seq)>1000000: val_msgs.append(f"⚠️ {seq_name}: Sequence very long (>{len(seq):,} bp)")
            if val_msgs:
                for msg in val_msgs: st.warning(msg)
                if any("Contains non-DNA characters" in m for m in val_msgs):
                    st.error("❌ Analysis stopped due to invalid sequence content."); st.session_state.analysis_status="Error"
            else:
                motif_results=[]
                with st.spinner("🧬 Analyzing motifs..."):
                    for i,seq in enumerate(st.session_state.seqs):
                        sequence_name=st.session_state.names[i] if i<len(st.session_state.names) else f"Sequence_{i+1}"
                        motifs=all_motifs_refactored(seq,sequence_name=sequence_name,nonoverlap=nonoverlap,report_hotspots=report_hotspots,calculate_conservation=calc_conservation)
                        validated_motifs=[ensure_subclass(m) for m in motifs if float(m.get('Normalized_Score',0) or 0)>=threshold]
                        motif_results.append(validated_motifs)
                st.session_state.results=motif_results
                summary=[]
                for i,motifs in enumerate(motif_results):
                    stats=get_basic_stats(st.session_state.seqs[i],motifs)
                    motif_types=Counter([m['Class'] if m['Class']!="Z-DNA" or m.get("Subclass")!="eGZ (Extruded-G)" else "eGZ (Extruded-G)" for m in motifs])
                    summary.append({"Sequence Name":st.session_state.names[i],"Length (bp)":stats['Length'],"GC %":stats['GC%'],"AT %":stats['AT%'],"A Count":stats['A'],"T Count":stats['T'],"G Count":stats['G'],"C Count":stats['C'],"Motif Count":len(motifs),"Motif Coverage (%)":stats.get("Motif Coverage %",0),"Motif Classes":", ".join(f"{k} ({v})" for k,v in motif_types.items())})
                st.session_state.summary_df=pd.DataFrame(summary)
                st.success("✅ Analysis complete! Results are available below and in the 'Analysis Results and Visualization' tab.")
                st.session_state.analysis_status="Complete"
        if st.session_state.get('summary_df') is not None:
            st.markdown("#### Analysis Summary"); st.dataframe(st.session_state.summary_df)
    st.markdown("---")

# --- Results, Download, Documentation Tabs: unchanged from your original for brevity ---
with tab_pages["Results"]:
    st.markdown('<h2>Analysis Results and Visualization</h2>',unsafe_allow_html=True)
    # ... (full logic unchanged, see your original file) ...
with tab_pages["Download"]:
    st.header("Export Data")
    # ... (full logic unchanged, see your original file) ...
with tab_pages["Documentation"]:
    st.header("Scientific Documentation & References")
    # ... (full logic unchanged, see your original file) ...

st.markdown("""
---
<div style='font-size: 1.05rem; color: #1e293b; margin-top: 30px; text-align: left; font-family:Montserrat,Arial;'>
<b>Developed by</b><br>
Dr. Venkata Rajesh Yella<br>
<a href='mailto:yvrajesh_bt@kluniversity.in'>yvrajesh_bt@kluniversity.in</a> |
<a href='https://github.com/VRYella' target='_blank'>GitHub: VRYella</a>
</div>
""",unsafe_allow_html=True)
