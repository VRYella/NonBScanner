# app.py - Streamlined NBDFinder Streamlit app
import io
import json
from collections import Counter
from typing import List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st
from Bio import Entrez, SeqIO

# Project helpers (assumed present in your repo)
from utils import parse_fasta, gc_content, reverse_complement, wrap
from motifs import get_basic_stats, visualization as viz
from motifs.enhanced_visualization import create_comprehensive_information_based_visualizations
from all_motifs_refactored import all_motifs_refactored

# Optional configuration module
try:
    from classification_config import MOTIF_LENGTH_LIMITS, SCORING_METHODS, get_motif_limits
    CONFIG_AVAILABLE = True
except Exception:
    CONFIG_AVAILABLE = False
    MOTIF_LENGTH_LIMITS = {}
    SCORING_METHODS = {}

# Optional export helpers
try:
    from nbdio.writers import export_to_bed, export_to_gff3, create_density_bedgraph
    NBDIO_AVAILABLE = True
except Exception:
    NBDIO_AVAILABLE = False

# ---------- APP CONFIG ----------
st.set_page_config(page_title="Non-B DNA Motif Finder", layout="wide", page_icon="🧬")
Entrez.email = "raazbiochem@gmail.com"

# ---------- SMALL HELPERS ----------
def ensure_subclass(motif: dict) -> dict:
    """Ensure motif has Subclass (or Subtype) set."""
    if not isinstance(motif, dict):
        return {"Subclass": "Other", "Motif": motif}
    if not motif.get("Subclass"):
        motif["Subclass"] = motif.get("Subtype", "Other")
    return motif

def parse_multi_fasta_text(content: str) -> Tuple[List[str], List[str]]:
    """Parse multi-FASTA plain text content into (seqs, names)."""
    seqs, names = [], []
    cur_seq, cur_name = "", ""
    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur_seq:
                seqs.append(parse_fasta(cur_seq))
                names.append(cur_name or f"Seq{len(seqs)}")
            cur_name = line.lstrip(">")
            cur_seq = ""
        else:
            cur_seq += line
    if cur_seq:
        seqs.append(parse_fasta(cur_seq))
        names.append(cur_name or f"Seq{len(seqs)}")
    return seqs, names

def show_sequence_preview(seqs: List[str], names: List[str], max_preview=2):
    """Show a compact preview for the first few sequences."""
    st.markdown("### 📊 Sequence Preview")
    for i, seq in enumerate(seqs[:max_preview]):
        stats = get_basic_stats(seq)
        st.markdown(
            f"**{names[i]}** ({len(seq):,} bp) | GC %: {stats['GC%']} | AT %: {stats['AT%']} | "
            f"A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}",
            unsafe_allow_html=True,
        )
        st.code(wrap(seq[:400]), language="fasta")
    if len(seqs) > max_preview:
        st.caption(f"...and {len(seqs)-max_preview} more.")

# Default UI strings/examples
EXAMPLE_FASTA = """>Example Sequence
ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC
"""
EXAMPLE_MULTI_FASTA = EXAMPLE_FASTA + "\n>Seq2\n" + EXAMPLE_FASTA.replace("Example Sequence", "Seq2")

# ---------- SESSION STATE DEFAULTS ----------
default_state = {
    "seqs": [],
    "names": [],
    "results": [],
    "summary_df": pd.DataFrame(),
    "analysis_status": "Ready",
}
for k, v in default_state.items():
    if k not in st.session_state:
        st.session_state[k] = v

# ---------- STYLES (light) ----------
st.markdown(
    """
    <style>
    .stApp { font-family: Montserrat, Arial, sans-serif; background: #f7fafd; }
    h1,h2,h3 { color: #1565c0; }
    .round-card { padding: 12px; border-radius:12px; background: #ffffff; box-shadow: 0 2px 6px rgba(0,0,0,0.04); }
    </style>
    """,
    unsafe_allow_html=True,
)

# ---------- PAGES & TABS ----------
PAGES = [
    "Home",
    "Upload & Analyze",
    "Results",
    "Download",
    "Documentation",
]
tabs = st.tabs(PAGES)
tab_home, tab_upload, tab_results, tab_download, tab_docs = tabs

# ---------- HOME ----------
with tab_home:
    st.header("Non-B DNA Motif Finder")
    c1, c2 = st.columns([1, 1])
    with c1:
        try:
            st.image("nbdcircle.JPG", use_column_width=True)
        except Exception:
            st.markdown(
                "<div class='round-card' style='text-align:center'>"
                "<h2>🧬 NBDFinder</h2><p>Non-B DNA motif detection and analysis.</p></div>",
                unsafe_allow_html=True,
            )
    with c2:
        st.markdown(
            """
            <div class='round-card'>
            <p><b>Detects 11 classes / 22+ subclasses of non-B DNA motifs.</b></p>
            <ul>
                <li>Upload multi-FASTA or paste sequences</li>
                <li>Run analysis (G4, Z-DNA, R-loops, Triplex, Slipped DNA, A-philic DNA, etc.)</li>
                <li>Download detailed tables, BED/GFF exports, and visualization suite</li>
            </ul>
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------- UPLOAD & ANALYZE ----------
with tab_upload:
    st.header("Sequence Upload and Motif Analysis")
    col_left, col_right = st.columns([1.1, 0.9], gap="large")

    # LEFT: Inputs & preview
    with col_left:
        st.subheader("Input Method")
        input_method = st.radio(
            "Choose input:",
            ["Upload FASTA", "Paste FASTA", "Example", "NCBI Fetch"],
            horizontal=True,
        )

        seqs, names = [], []
        if input_method == "Upload FASTA":
            fasta_file = st.file_uploader("Drop FASTA / multi-FASTA", type=["fa", "fasta", "txt"])
            if fasta_file:
                content = fasta_file.read().decode("utf-8")
                seqs, names = parse_multi_fasta_text(content)
                st.success(f"Loaded {len(seqs)} sequence(s).")

        elif input_method == "Paste FASTA":
            text = st.text_area("Paste single or multi-FASTA", height=140)
            if text:
                seqs, names = parse_multi_fasta_text(text)
                st.success(f"Pasted {len(seqs)} sequence(s).")

        elif input_method == "Example":
            ex_type = st.radio("Example type", ["Single", "Multi"], horizontal=True)
            if st.button("Load Example"):
                if ex_type == "Single":
                    seqs, names = [parse_fasta(EXAMPLE_FASTA)], ["Example Sequence"]
                else:
                    seqs, names = parse_multi_fasta_text(EXAMPLE_MULTI_FASTA)
                st.success("Example loaded.")

        else:  # NCBI Fetch
            db = st.selectbox("NCBI DB", ["nucleotide", "gene"])
            query = st.text_input("Accession / query")
            rettype = st.radio("Return format", ["fasta", "gb"], horizontal=True)
            retmax = st.number_input("Max records", min_value=1, max_value=20, value=3)
            motif_examples = {
                "G-quadruplex": "NR_003287.2",
                "Z-DNA": "NM_001126112.2",
                "R-loop": "NR_024540.1",
            }
            with st.expander("Motif example queries"):
                for k, v in motif_examples.items():
                    st.write(f"**{k}**: `{v}`")
            if st.button("Fetch from NCBI") and query:
                with st.spinner("Fetching..."):
                    try:
                        handle = Entrez.efetch(db=db, id=query, rettype=rettype, retmode="text")
                        records = list(SeqIO.parse(handle, rettype))
                        handle.close()
                        seqs = [str(rec.seq).upper().replace("U", "T") for rec in records]
                        names = [rec.id for rec in records]
                        st.success(f"Fetched {len(seqs)} record(s).")
                    except Exception as e:
                        st.error(f"NCBI fetch error: {e}")

        # Persist into session state if we got sequences
        if seqs:
            st.session_state.seqs = seqs
            st.session_state.names = names
            st.session_state.results = []
            show_sequence_preview(seqs, names)

        # If session already has sequences (from previous action)
        elif st.session_state.seqs:
            show_sequence_preview(st.session_state.seqs, st.session_state.names)

    # RIGHT: Run analysis control + short config
    with col_right:
        st.subheader("Run Analysis")
        st.markdown(
            f"""
            **Configuration summary:**  
            - Motif classes: {len(MOTIF_LENGTH_LIMITS) if CONFIG_AVAILABLE else 'default set'}  
            - Scoring: normalized (0–1) + raw scores  
            - Threshold: include all detected motifs (internal threshold = 0.0)
            """
        )

        # Run button uses safe defaults (no UI advanced controls shown)
        if st.button("🔬 Run Motif Analysis", use_container_width=True):
            if not st.session_state.seqs:
                st.warning("No sequences available. Upload or paste sequences first.")
            else:
                st.session_state.analysis_status = "Running"
                validation_messages = []
                # basic sequence validation
                for i, seq in enumerate(st.session_state.seqs):
                    name = st.session_state.names[i] if i < len(st.session_state.names) else f"Sequence_{i+1}"
                    valid_chars = set("ATCGN")
                    seq_chars = set(seq.upper())
                    if not seq_chars.issubset(valid_chars):
                        invalid = seq_chars - valid_chars
                        validation_messages.append(f"{name}: non DNA characters {invalid}")
                    if len(seq) < 10:
                        validation_messages.append(f"{name}: too short (<10 bp)")

                if validation_messages:
                    for m in validation_messages:
                        st.warning(m)
                    st.session_state.analysis_status = "Error"
                else:
                    motif_results = []
                    with st.spinner("Analyzing motifs..."):
                        for i, seq in enumerate(st.session_state.seqs):
                            name = st.session_state.names[i] if i < len(st.session_state.names) else f"Sequence_{i+1}"
                            motifs = all_motifs_refactored(
                                seq,
                                sequence_name=name,
                                nonoverlap=False,
                                report_hotspots=True,
                                calculate_conservation=False,
                            )
                            # Normalize and ensure subclass
                            validated = []
                            for m in motifs:
                                try:
                                    m_score = float(m.get("Normalized_Score", 0))
                                except Exception:
                                    m_score = 0.0
                                if m_score >= 0.0:
                                    validated.append(ensure_subclass(m))
                            motif_results.append(validated)
                    st.session_state.results = motif_results
                    # Build summary dataframe
                    rows = []
                    for i, motifs in enumerate(motif_results):
                        stats = get_basic_stats(st.session_state.seqs[i], motifs)
                        cls_counts = Counter([m.get("Class", "Unknown") for m in motifs])
                        rows.append(
                            {
                                "Sequence Name": st.session_state.names[i],
                                "Length (bp)": stats["Length"],
                                "GC %": stats["GC%"],
                                "Motif Count": len(motifs),
                                "Motif Coverage (%)": stats.get("Motif Coverage %", 0),
                                "Motif Classes": ", ".join(f"{k} ({v})" for k, v in cls_counts.items()),
                            }
                        )
                    st.session_state.summary_df = pd.DataFrame(rows)
                    st.session_state.analysis_status = "Complete"
                    st.success("Analysis complete. See Results tab.")

# ---------- RESULTS ----------
with tab_results:
    st.header("Analysis Results and Visualization")
    if not st.session_state.results:
        st.info("No results. Run analysis from 'Upload & Analyze'.")
    else:
        st.subheader("Summary")
        st.dataframe(st.session_state.summary_df, use_container_width=True)

        # Select sequence to explore
        idx = 0
        if len(st.session_state.seqs) > 1:
            idx = st.selectbox(
                "Choose sequence", options=list(range(len(st.session_state.seqs))),
                format_func=lambda i: st.session_state.names[i]
            )

        motifs = st.session_state.results[idx]
        seq = st.session_state.seqs[idx]
        seq_name = st.session_state.names[idx]
        seq_len = len(seq)

        if not motifs:
            st.warning("No motifs found for this sequence.")
        else:
            # DataFrame of motifs
            df = pd.DataFrame(motifs)

            # Quick metrics card
            filtered = [m for m in motifs if m.get("Class") not in ("Hybrid", "Non-B DNA Clusters", "Non-B_DNA_Cluster")]
            excluded = [m for m in motifs if m.get("Class") in ("Hybrid", "Non-B DNA Clusters", "Non-B_DNA_Cluster")]
            coverage = get_basic_stats(seq, motifs).get("Motif Coverage %", 0)
            non_b_density = (len(filtered) / seq_len * 1000) if seq_len else 0

            st.markdown(
                f"""
                <div class='round-card'>
                    <b>Sequence:</b> {seq_name} &nbsp; | &nbsp; <b>Length:</b> {seq_len:,} bp
                    <br><b>Motifs:</b> {len(motifs)} &nbsp; | &nbsp; <b>Coverage:</b> {coverage:.2f}% &nbsp; | &nbsp; <b>Density:</b> {non_b_density:.2f} motifs/kb
                </div>""",
                unsafe_allow_html=True,
            )

            # Score analysis (histograms) if available
            if any("Normalized_Score" in m or "Score" in m for m in motifs):
                st.subheader("Score Distributions")
                fig, axes = plt.subplots(1, 2, figsize=(10, 4))
                actual = [float(m.get("Actual_Score", m.get("Score", 0)) or 0) for m in motifs]
                normed = [float(m.get("Normalized_Score", 0) or 0) for m in motifs]
                if any(actual):
                    axes[0].hist(actual, bins=20)
                    axes[0].set_title("Actual Scores")
                else:
                    axes[0].text(0.5, 0.5, "No actual scores", ha="center")
                if any(normed):
                    axes[1].hist(normed, bins=20)
                    axes[1].set_title("Normalized Scores")
                else:
                    axes[1].text(0.5, 0.5, "No normalized scores", ha="center")
                plt.tight_layout()
                st.pyplot(fig)
                plt.close(fig)

            # Motif class distribution
            st.subheader("Motif Class Distribution")
            class_counts = df["Class"].value_counts()
            fig, ax = plt.subplots(figsize=(8, max(3, 0.25 * len(class_counts))))
            ax.barh(class_counts.index, class_counts.values)
            ax.set_xlabel("Count")
            st.pyplot(fig)
            plt.close(fig)

            # Detailed motif table: allow column selection
            st.subheader("Detailed motif table")
            available_cols = df.columns.tolist()
            default_cols = [c for c in ["Class", "Subclass", "Start", "End", "Length", "Normalized_Score", "Actual_Score", "GC_Content"] if c in available_cols]
            sel = st.multiselect("Columns to display", available_cols, default=default_cols)
            disp_df = df[sel].copy() if sel else df.copy()
            # Clean up column names for display
            disp_df.columns = [c.replace("_", " ") for c in disp_df.columns]
            st.dataframe(disp_df, use_container_width=True, height=360)

            # Comprehensive visualization generation (best-effort)
            st.subheader("Comprehensive visualizations")
            with st.spinner("Generating visualizations..."):
                try:
                    static_plots, interactive_plots, _ = create_comprehensive_information_based_visualizations(df, seq_len, seq_name)
                    # show a few categories if available
                    for key, fig in static_plots.items():
                        st.markdown(f"**{key.replace('_', ' ').title()}**")
                        st.pyplot(fig)
                        plt.close(fig)
                    if interactive_plots:
                        for name, fig in interactive_plots.items():
                            st.plotly_chart(fig, use_container_width=True)
                except Exception as e:
                    st.error(f"Visualization error: {e}")
                    # Fallback simple plot: motif position map
                    st.info("Using fallback motif map")
                    fig, ax = plt.subplots(figsize=(12, 2 + 0.2 * len(df)))
                    for i, (_, r) in enumerate(df.iterrows()):
                        start = r.get("Start", 0)
                        end = r.get("End", start + 1)
                        ax.hlines(i, start, end, lw=4)
                    ax.set_xlabel("Position (bp)")
                    ax.set_ylabel("Motif index")
                    st.pyplot(fig)
                    plt.close(fig)

# ---------- DOWNLOAD ----------
with tab_download:
    st.header("Export Data")
    if not st.session_state.results:
        st.info("No results to export.")
    else:
        # Flatten results
        rows = []
        for i, motifs in enumerate(st.session_state.results):
            for m in motifs:
                r = m.copy()
                r.setdefault("Sequence Name", st.session_state.names[i])
                # Normalize class name display
                if r.get("Class") == "Z-DNA" and r.get("Subclass") == "eGZ (Extruded-G)":
                    r["Class"] = "eGZ (Extruded-G)"
                rows.append(r)
        df_all = pd.DataFrame(rows).loc[:, ~pd.DataFrame(rows).columns.duplicated()] if rows else pd.DataFrame()

        st.subheader("Export preview")
        st.dataframe(df_all.head(10), use_container_width=True)

        col1, col2, col3 = st.columns(3)
        with col1:
            csv_bytes = df_all.to_csv(index=False).encode("utf-8")
            st.download_button("Download CSV", data=csv_bytes, file_name="nbdfinder_results.csv", mime="text/csv")
        with col2:
            out = io.BytesIO()
            with pd.ExcelWriter(out, engine="xlsxwriter") as writer:
                df_all.to_excel(writer, index=False, sheet_name="Motifs")
            out.seek(0)
            st.download_button("Download Excel", data=out, file_name="nbdfinder_results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        with col3:
            if CONFIG_AVAILABLE:
                conf = {
                    "motif_length_limits": MOTIF_LENGTH_LIMITS,
                    "scoring_methods": SCORING_METHODS,
                }
                st.download_button("Download Config JSON", data=json.dumps(conf, indent=2), file_name="nbdfinder_config.json", mime="application/json")

        # Genome browser exports (optional)
        if NBDIO_AVAILABLE:
            st.subheader("Genome browser formats")
            col4, col5, col6 = st.columns(3)
            all_motifs_flat = [m for seq_m in st.session_state.results for m in seq_m]
            with col4:
                try:
                    bed = export_to_bed(all_motifs_flat, sequence_name="NBDFinder_Analysis", score_type="normalized")
                    st.download_button("Download BED", data=bed, file_name="nbdfinder.bed", mime="text/plain")
                except Exception as e:
                    st.error(f"BED export error: {e}")
            with col5:
                try:
                    gff3 = export_to_gff3(all_motifs_flat, sequence_name="NBDFinder_Analysis")
                    st.download_button("Download GFF3", data=gff3, file_name="nbdfinder.gff3", mime="text/plain")
                except Exception as e:
                    st.error(f"GFF3 export error: {e}")
            with col6:
                try:
                    seq_len_est = max((m.get("End", 0) for m in all_motifs_flat), default=1000)
                    density = create_density_bedgraph([m for m in all_motifs_flat if m.get("Class") not in ("Hybrid", "Non-B_DNA_Cluster")], seq_len_est, sequence_name="NBDFinder_Analysis")
                    st.download_button("Download density bedGraph", data=density, file_name="nbdfinder_density.bedgraph", mime="text/plain")
                except Exception as e:
                    st.error(f"Density export error: {e}")
        else:
            st.info("Genome browser export support not installed (nbdio).")

# ---------- DOCUMENTATION ----------
with tab_docs:
    st.header("Documentation & References")
    st.markdown(
        """
        - Detects canonical and non-canonical DNA motif types (G4, Z-DNA, R-loop, Triplex, Slipped DNA, A-philic DNA, etc.)
        - Visualizations: coverage, density, class distributions, motif maps, interactive plots (if available)
        - API: possibility for REST API (see repo docs)
        """
    )

# Footer / author
st.markdown(
    """--- 
    <div style='font-size:0.95rem'>
    Developed by Dr. Venkata Rajesh Yella — <a href='mailto:yvrajesh_bt@kluniversity.in'>yvrajesh_bt@kluniversity.in</a>
    </div>
    """,
    unsafe_allow_html=True,
)
