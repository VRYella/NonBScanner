"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                         NBDSCANNER WEB APPLICATION                            ║
║                    Non-B DNA Motif Detection System                          ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: app.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    Streamlit web application for comprehensive Non-B DNA motif detection.
    Provides interactive interface for sequence analysis and visualization.

FEATURES:
┌─────────────────────────────────────────────────────────────────────────────┐
│  - Multi-FASTA support             - Real-time analysis progress            │
│  - 11 motif classes detection      - Interactive visualizations             │
│  - 22+ subclass analysis           - Export to CSV/BED/JSON                 │
│  - NCBI sequence fetch             - Publication-quality plots              │
└─────────────────────────────────────────────────────────────────────────────┘

ARCHITECTURE:
    Input → Detection → Scoring → Overlap Resolution → Visualization → Export
"""

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import io
import numpy as np
from collections import Counter
# Import consolidated NBDScanner modules
from utilities import (
    canonicalize_motif, parse_fasta, gc_content, reverse_complement, wrap,
    get_basic_stats, export_to_bed, export_to_csv, export_to_json,
    validate_sequence, quality_check_motifs
)
from scanner import (
    analyze_sequence, analyze_multiple_sequences,
    get_motif_classification_info, export_results_to_dataframe
)
from visualizations import (
    plot_motif_distribution, plot_coverage_map,
    plot_length_distribution, plot_nested_pie_chart, save_all_plots,
    MOTIF_CLASS_COLORS
)

# Try to import Entrez for demo functionality
try:
    from Bio import Entrez, SeqIO
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False

# ---------- PAGE CONFIG ----------
st.set_page_config(
    page_title="NBDScanner - Non-B DNA Motif Finder",
    layout="wide",
    page_icon="🧬",
    menu_items={'About': "NBDScanner | Developed by Dr. Venkata Rajesh Yella"}
)

# Get motif classification info
CLASSIFICATION_INFO = get_motif_classification_info()

# ---------- PATCH: Ensure every motif has Subclass ----------
def ensure_subclass(motif):
    """Guarantee every motif has a string 'Subclass'"""
    if isinstance(motif, dict):
        if 'Subclass' not in motif or motif['Subclass'] is None:
            motif['Subclass'] = motif.get('Subtype', 'Other')
        return motif
    else:
        # Handle non-dict motifs gracefully
        return {'Subclass': 'Other', 'Motif': motif}


# ---------- ENHANCED PROFESSIONAL CSS FOR RESEARCH-QUALITY UI ----------
st.markdown("""
    <style>
    /* Import Google Fonts for professional scientific typography */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=IBM+Plex+Sans:wght@300;400;500;600;700&family=Source+Sans+Pro:wght@400;600;700&display=swap');
    
    /* Main app background with elegant scientific gradient */
    body, [data-testid="stAppViewContainer"], .main {
        background: linear-gradient(135deg, #f5f9ff 0%, #e8f4fd 50%, #f0f7ff 100%) !important;
        font-family: 'Inter', 'IBM Plex Sans', 'Segoe UI', system-ui, -apple-system, sans-serif !important;
    }
    
    /* Tabs: Premium scientific design with advanced styling */
    .stTabs [data-baseweb="tab-list"] {
        width: 98vw !important;
        justify-content: stretch !important;
        border-bottom: 3px solid #0d47a1;
        background: linear-gradient(90deg, #ffffff 0%, #e3f2fd 50%, #f8fbff 100%) !important;
        box-shadow: 0 6px 16px rgba(13, 71, 161, 0.12);
        margin-bottom: 1.5em;
        border-radius: 12px 12px 0 0;
    }
    .stTabs [data-baseweb="tab"] {
        font-size: 1.1rem !important;
        font-weight: 600 !important;
        flex: 1 1 0%;
        min-width: 0 !important;
        padding: 16px 12px !important;
        text-align: center;
        color: #455a64 !important;
        background: transparent !important;
        border-right: 1px solid rgba(0,0,0,0.05) !important;
        letter-spacing: 0.04em;
        transition: all 0.35s cubic-bezier(0.4, 0, 0.2, 1);
        position: relative;
    }
    .stTabs [data-baseweb="tab"]::after {
        content: '';
        position: absolute;
        bottom: 0;
        left: 50%;
        width: 0%;
        height: 4px;
        background: linear-gradient(90deg, #0d47a1 0%, #1976d2 100%);
        transition: all 0.35s cubic-bezier(0.4, 0, 0.2, 1);
        transform: translateX(-50%);
        border-radius: 4px 4px 0 0;
    }
    .stTabs [data-baseweb="tab"]:hover {
        background: linear-gradient(180deg, rgba(13, 71, 161, 0.05) 0%, rgba(25, 118, 210, 0.08) 100%) !important;
        color: #0d47a1 !important;
        transform: translateY(-2px);
    }
    .stTabs [data-baseweb="tab"]:hover::after {
        width: 60%;
    }
    .stTabs [aria-selected="true"] {
        color: #0d47a1 !important;
        background: linear-gradient(180deg, #ffffff 0%, #e3f2fd 100%) !important;
        box-shadow: 0 -4px 12px rgba(13, 71, 161, 0.15), inset 0 -3px 0 0 #0d47a1;
        font-weight: 700 !important;
        transform: translateY(-2px);
    }
    .stTabs [aria-selected="true"]::after {
        width: 100%;
        height: 5px;
    }
    .stTabs [data-baseweb="tab"]:last-child {
        border-right: none !important;
    }
    
    /* Headings: Enhanced scientific typography hierarchy */
    h1, h2, h3, h4 {
        font-family: 'IBM Plex Sans', 'Inter', 'Segoe UI', system-ui, sans-serif !important;
        color: #0d47a1 !important;
        font-weight: 700 !important;
        letter-spacing: -0.025em;
        margin-top: 1.3em;
        margin-bottom: 0.8em;
        text-shadow: 0 1px 2px rgba(0, 0, 0, 0.05);
    }
    h1 { 
        font-size: 2.4rem !important; 
        font-weight: 800 !important;
        background: linear-gradient(135deg, #0d47a1 0%, #1976d2 50%, #42a5f5 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
        letter-spacing: -0.03em;
    }
    h2 { 
        font-size: 1.7rem !important; 
        color: #0d47a1 !important; 
        font-weight: 700 !important;
        border-bottom: 3px solid #e3f2fd;
        padding-bottom: 0.5rem;
    }
    h3 { 
        font-size: 1.3rem !important; 
        color: #1565c0 !important; 
        font-weight: 700 !important;
    }
    h4 { 
        font-size: 1.15rem !important; 
        color: #1976d2 !important; 
        font-weight: 600 !important;
    }
    
    /* Body text: Scientific readability with optimal contrast */
    .stMarkdown, .markdown-text-container, .stText, p, span, label {
        font-size: 1.02rem !important;
        font-family: 'Source Sans Pro', 'Inter', 'Segoe UI', system-ui, sans-serif !important;
        line-height: 1.75 !important;
        color: #1e293b !important;
        margin-bottom: 0.75em !important;
        font-weight: 400;
    }
    
    /* Input fields: Premium scientific design with elegant borders */
    input, .stTextInput>div>div>input, .stSelectbox>div>div>div, 
    .stMultiSelect>div>div>div, .stRadio>div>div>label>div {
        font-size: 1.0rem !important;
        font-family: 'Inter', 'IBM Plex Sans', 'Segoe UI', system-ui, sans-serif !important;
        border-radius: 10px !important;
        border: 2px solid #e1e8ed !important;
        background: #ffffff !important;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.02);
    }
    input:focus, .stTextInput>div>div>input:focus {
        border-color: #1976d2 !important;
        box-shadow: 0 0 0 4px rgba(25, 118, 210, 0.12), 0 4px 8px rgba(25, 118, 210, 0.08) !important;
        background: #fefeff !important;
        transform: translateY(-1px);
    }
    
    /* Buttons: Advanced gradient design with scientific elegance */
    .stButton>button {
        font-size: 1.05rem !important;
        font-family: 'Inter', 'IBM Plex Sans', 'Segoe UI', system-ui, sans-serif !important;
        padding: 0.75em 1.8em !important;
        background: linear-gradient(135deg, #0d47a1 0%, #1976d2 50%, #42a5f5 100%) !important;
        color: #fff !important;
        border-radius: 12px !important;
        border: none !important;
        font-weight: 600 !important;
        box-shadow: 0 6px 18px rgba(13, 71, 161, 0.35), 0 2px 6px rgba(25, 118, 210, 0.2);
        transition: all 0.35s cubic-bezier(0.4, 0, 0.2, 1);
        letter-spacing: 0.03em;
        position: relative;
        overflow: hidden;
    }
    .stButton>button::before {
        content: '';
        position: absolute;
        top: 0;
        left: -100%;
        width: 100%;
        height: 100%;
        background: linear-gradient(90deg, transparent, rgba(255,255,255,0.3), transparent);
        transition: left 0.5s;
    }
    .stButton>button:hover::before {
        left: 100%;
    }
    .stButton>button:hover {
        background: linear-gradient(135deg, #0a3d91 0%, #1565c0 50%, #1976d2 100%) !important;
        box-shadow: 0 8px 24px rgba(13, 71, 161, 0.45), 0 4px 12px rgba(25, 118, 210, 0.3);
        transform: translateY(-3px) scale(1.02);
    }
    .stButton>button:active {
        transform: translateY(-1px) scale(0.98);
        box-shadow: 0 4px 12px rgba(13, 71, 161, 0.35), 0 2px 6px rgba(25, 118, 210, 0.2);
    }
    
    /* DataFrames: Scientific table design with elegant borders */
    .stDataFrame, .stTable {
        font-size: 0.96rem !important;
        font-family: 'IBM Plex Sans', 'Inter', 'Segoe UI', system-ui, sans-serif !important;
        line-height: 1.65 !important;
        border-radius: 12px !important;
        overflow: hidden;
        box-shadow: 0 4px 12px rgba(13, 71, 161, 0.08), 0 2px 4px rgba(0, 0, 0, 0.03) !important;
        border: 1px solid #e3f2fd !important;
    }
    .stDataFrame thead tr th {
        background: linear-gradient(135deg, #0d47a1 0%, #1976d2 100%) !important;
        color: white !important;
        font-weight: 700 !important;
        padding: 1rem !important;
        font-size: 0.98rem !important;
        letter-spacing: 0.02em;
        text-transform: uppercase;
        border-bottom: 3px solid #0a3d91 !important;
    }
    .stDataFrame tbody tr {
        transition: all 0.25s ease;
        border-bottom: 1px solid #e3f2fd !important;
    }
    .stDataFrame tbody tr:hover {
        background: linear-gradient(90deg, #f0f7ff 0%, #e3f2fd 100%) !important;
        transform: scale(1.005);
        box-shadow: 0 2px 8px rgba(25, 118, 210, 0.1);
    }
    .stDataFrame tbody tr td {
        padding: 0.85rem !important;
        color: #1e293b !important;
    }
    
    /* Tab content spacing */
    .stTabs [data-baseweb="tab-panel"] {
        padding-top: 2.5rem !important;
        padding-left: 1.5rem !important;
        padding-right: 1.5rem !important;
        padding-bottom: 2rem !important;
    }
    
    /* Analysis summary cards with modern shadows */
    .analysis-summary-card {
        margin: 1.5rem 0 !important;
        padding: 2rem !important;
        background: white !important;
        border-radius: 12px !important;
        box-shadow: 0 4px 16px rgba(0,0,0,0.08) !important;
        border: 1px solid rgba(25, 118, 210, 0.08) !important;
    }
    
    /* Select boxes and multiselect: Premium dropdown design */
    .stSelectbox > div > div > div, .stMultiSelect > div > div > div {
        min-height: 3.5rem !important;
        padding: 1rem 1.2rem !important;
        line-height: 1.65 !important;
        border-radius: 12px !important;
        background: linear-gradient(135deg, #ffffff 0%, #f8fbff 100%) !important;
        border: 2px solid #e1e8ed !important;
        box-shadow: 0 2px 6px rgba(0, 0, 0, 0.03);
        transition: all 0.3s ease;
    }
    .stSelectbox > div > div > div:hover, .stMultiSelect > div > div > div:hover {
        border-color: #64b5f6 !important;
        box-shadow: 0 4px 12px rgba(100, 181, 246, 0.12);
        transform: translateY(-1px);
    }
    
    /* Enhanced input field spacing with elegant styling */
    .stTextInput > div > div > input, .stTextArea textarea {
        padding: 1rem 1.3rem !important;
        min-height: 3.5rem !important;
        line-height: 1.6 !important;
        border-radius: 12px !important;
        font-size: 1.02rem !important;
        background: linear-gradient(135deg, #ffffff 0%, #f8fbff 100%) !important;
    }
    .stTextArea textarea {
        min-height: 120px !important;
    }
    
    /* Number input styling with scientific precision */
    .stNumberInput > div > div > input {
        padding: 1rem 1.3rem !important;
        min-height: 3.5rem !important;
        border-radius: 12px !important;
        font-size: 1.02rem !important;
        font-weight: 600;
        background: linear-gradient(135deg, #ffffff 0%, #f8fbff 100%) !important;
    }
    
    /* Radio buttons: Enhanced with elegant hover states and borders */
    .stRadio > div {
        gap: 0.4rem !important;
    }
    .stRadio > div > div > label {
        margin-bottom: 0.7rem !important;
        padding: 0.65rem 1rem !important;
        border-radius: 10px !important;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        border: 2px solid transparent;
        background: linear-gradient(135deg, #ffffff 0%, #f8fbff 100%);
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.03);
        cursor: pointer;
    }
    .stRadio > div > div > label:hover {
        background: linear-gradient(135deg, #e3f2fd 0%, #f0f7ff 100%) !important;
        border-color: #bbdefb !important;
        transform: translateX(4px);
        box-shadow: 0 4px 8px rgba(25, 118, 210, 0.08);
    }
    .stRadio > div > div > label[data-checked="true"] {
        background: linear-gradient(135deg, #e3f2fd 0%, #bbdefb 100%) !important;
        border-color: #1976d2 !important;
        font-weight: 600 !important;
        box-shadow: 0 4px 12px rgba(25, 118, 210, 0.15);
    }
    /* Radio button circles */
    .stRadio > div > div > label > div:first-child {
        border-width: 2px !important;
        border-color: #90caf9 !important;
        width: 22px !important;
        height: 22px !important;
    }
    .stRadio > div > div > label[data-checked="true"] > div:first-child {
        border-color: #0d47a1 !important;
        background-color: #0d47a1 !important;
    }
    
    /* Number input styling */
    .stNumberInput > div > div > input {
        padding: 0.85rem 1.1rem !important;
        min-height: 3rem !important;
        border-radius: 8px !important;
    }
    
    /* Info boxes and alerts: Scientific notification design */
    .stAlert {
        border-radius: 12px !important;
        border-left: 5px solid #1976d2 !important;
        box-shadow: 0 4px 12px rgba(0,0,0,0.08) !important;
        background: linear-gradient(90deg, #e3f2fd 0%, #f0f7ff 100%) !important;
        padding: 1.2rem !important;
        font-weight: 500;
    }
    .stAlert[data-baseweb="notification"] {
        border-left-color: #0d47a1 !important;
    }
    .stSuccess {
        border-left-color: #2e7d32 !important;
        background: linear-gradient(90deg, #e8f5e9 0%, #f1f8f4 100%) !important;
    }
    .stWarning {
        border-left-color: #f57c00 !important;
        background: linear-gradient(90deg, #fff3e0 0%, #fff8f1 100%) !important;
    }
    .stError {
        border-left-color: #c62828 !important;
        background: linear-gradient(90deg, #ffebee 0%, #fff5f5 100%) !important;
    }
    
    /* Progress bars: Elegant animated design */
    .stProgress > div > div {
        background: linear-gradient(90deg, #0d47a1 0%, #1976d2 50%, #42a5f5 100%) !important;
        border-radius: 12px !important;
        box-shadow: 0 2px 8px rgba(25, 118, 210, 0.3);
        animation: shimmer 2s infinite;
    }
    @keyframes shimmer {
        0% { background-position: -1000px 0; }
        100% { background-position: 1000px 0; }
    }
    .stProgress > div {
        border-radius: 12px !important;
        background: #e3f2fd !important;
        box-shadow: inset 0 2px 4px rgba(0, 0, 0, 0.05);
    }
    
    /* File uploader: Premium drag-and-drop design */
    [data-testid="stFileUploader"] {
        border: 3px dashed #64b5f6 !important;
        border-radius: 16px !important;
        background: linear-gradient(135deg, rgba(227, 242, 253, 0.3) 0%, rgba(187, 222, 251, 0.2) 100%) !important;
        padding: 2.5rem !important;
        transition: all 0.35s cubic-bezier(0.4, 0, 0.2, 1);
        box-shadow: 0 4px 12px rgba(100, 181, 246, 0.1);
    }
    [data-testid="stFileUploader"]:hover {
        border-color: #0d47a1 !important;
        background: linear-gradient(135deg, rgba(227, 242, 253, 0.5) 0%, rgba(187, 222, 251, 0.4) 100%) !important;
        transform: scale(1.01);
        box-shadow: 0 6px 20px rgba(13, 71, 161, 0.15);
    }
    [data-testid="stFileUploader"] section {
        border: none !important;
    }
    [data-testid="stFileUploader"] button {
        background: linear-gradient(135deg, #1976d2 0%, #42a5f5 100%) !important;
        color: white !important;
        border: none !important;
        border-radius: 10px !important;
        padding: 0.6rem 1.5rem !important;
        font-weight: 600 !important;
        transition: all 0.3s ease;
    }
    [data-testid="stFileUploader"] button:hover {
        background: linear-gradient(135deg, #0d47a1 0%, #1976d2 100%) !important;
        transform: translateY(-2px);
        box-shadow: 0 6px 16px rgba(25, 118, 210, 0.3);
    }
    
    /* Checkboxes: Modern toggle-style design */
    .stCheckbox {
        padding: 0.5rem 0 !important;
    }
    .stCheckbox > label {
        padding: 0.6rem 0.8rem !important;
        border-radius: 10px !important;
        transition: all 0.3s ease;
        cursor: pointer;
        background: linear-gradient(135deg, #ffffff 0%, #f8fbff 100%);
        border: 2px solid transparent;
    }
    .stCheckbox > label:hover {
        background: linear-gradient(135deg, #e3f2fd 0%, #f0f7ff 100%) !important;
        border-color: #bbdefb !important;
        transform: translateX(3px);
    }
    .stCheckbox > label > div:first-child {
        border-width: 2px !important;
        border-color: #64b5f6 !important;
        border-radius: 6px !important;
        width: 22px !important;
        height: 22px !important;
        transition: all 0.3s ease;
    }
    .stCheckbox > label > div:first-child[data-checked="true"] {
        background-color: #0d47a1 !important;
        border-color: #0d47a1 !important;
    }
    
    /* Expander styling: Scientific accordion design */
    .streamlit-expanderHeader {
        border-radius: 12px !important;
        background: linear-gradient(135deg, #e3f2fd 0%, #f0f7ff 100%) !important;
        font-weight: 700 !important;
        padding: 1rem 1.5rem !important;
        border: 2px solid #bbdefb !important;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        box-shadow: 0 2px 6px rgba(25, 118, 210, 0.08);
    }
    .streamlit-expanderHeader:hover {
        background: linear-gradient(135deg, #bbdefb 0%, #e3f2fd 100%) !important;
        border-color: #1976d2 !important;
        transform: translateX(4px);
        box-shadow: 0 4px 12px rgba(25, 118, 210, 0.15);
    }
    .streamlit-expanderContent {
        border-radius: 0 0 12px 12px !important;
        background: #fefeff !important;
        border: 2px solid #e3f2fd !important;
        border-top: none !important;
        padding: 1.5rem !important;
    }
    
    /* Metric cards: Advanced scientific metrics */
    [data-testid="stMetric"] {
        background: linear-gradient(135deg, #ffffff 0%, #f0f7ff 100%) !important;
        padding: 1.5rem !important;
        border-radius: 14px !important;
        border: 2px solid #e3f2fd !important;
        box-shadow: 0 4px 12px rgba(13, 71, 161, 0.08) !important;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
    }
    [data-testid="stMetric"]:hover {
        transform: translateY(-4px) scale(1.02);
        box-shadow: 0 8px 24px rgba(13, 71, 161, 0.15) !important;
        border-color: #64b5f6 !important;
    }
    [data-testid="stMetricValue"] {
        font-size: 2.2rem !important;
        font-weight: 800 !important;
        background: linear-gradient(135deg, #0d47a1 0%, #1976d2 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
    }
    [data-testid="stMetricLabel"] {
        font-size: 0.95rem !important;
        font-weight: 600 !important;
        color: #455a64 !important;
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }
    
    /* Spinner */
    .stSpinner > div {
        border-top-color: #1976d2 !important;
    }
    
    /* Code blocks: Scientific monospace design */
    .stCodeBlock {
        border-radius: 12px !important;
        border: 2px solid #e3f2fd !important;
        background: linear-gradient(135deg, #f8fbff 0%, #f0f7ff 100%) !important;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.03);
        padding: 1rem !important;
    }
    code {
        font-family: 'JetBrains Mono', 'Fira Code', 'Consolas', monospace !important;
        font-size: 0.92rem !important;
        color: #0d47a1 !important;
        background: rgba(227, 242, 253, 0.4) !important;
        padding: 0.2rem 0.5rem !important;
        border-radius: 6px !important;
    }
    
    /* Sidebar styling for navigation */
    [data-testid="stSidebar"] {
        background: linear-gradient(180deg, #e3f2fd 0%, #f0f7ff 100%) !important;
        border-right: 2px solid #bbdefb !important;
    }
    [data-testid="stSidebar"] .stMarkdown {
        color: #0d47a1 !important;
    }
    
    /* Scrollbar styling: Elegant scientific design */
    ::-webkit-scrollbar {
        width: 12px;
        height: 12px;
    }
    ::-webkit-scrollbar-track {
        background: linear-gradient(135deg, #f1f5f9 0%, #e3f2fd 100%);
        border-radius: 12px;
        border: 1px solid #e3f2fd;
    }
    ::-webkit-scrollbar-thumb {
        background: linear-gradient(135deg, #0d47a1 0%, #1976d2 50%, #42a5f5 100%);
        border-radius: 12px;
        border: 2px solid #f8fbff;
        box-shadow: 0 2px 6px rgba(13, 71, 161, 0.2);
    }
    ::-webkit-scrollbar-thumb:hover {
        background: linear-gradient(135deg, #0a3d91 0%, #1565c0 50%, #1976d2 100%);
        box-shadow: 0 3px 10px rgba(13, 71, 161, 0.3);
    }
    
    /* Loading spinner */
    .stSpinner > div {
        border-top-color: #0d47a1 !important;
        border-right-color: #1976d2 !important;
        border-bottom-color: #42a5f5 !important;
    }
    
    /* Download buttons enhanced */
    [data-testid="stDownloadButton"] button {
        background: linear-gradient(135deg, #0d47a1 0%, #1976d2 100%) !important;
        color: white !important;
        border-radius: 12px !important;
        padding: 0.7rem 1.6rem !important;
        font-weight: 600 !important;
        box-shadow: 0 4px 12px rgba(13, 71, 161, 0.25);
        transition: all 0.3s ease;
    }
    [data-testid="stDownloadButton"] button:hover {
        background: linear-gradient(135deg, #0a3d91 0%, #1565c0 100%) !important;
        transform: translateY(-2px);
        box-shadow: 0 6px 18px rgba(13, 71, 161, 0.35);
    }
    
    /* Caption text styling */
    .caption, [data-testid="stCaptionContainer"] {
        color: #607d8b !important;
        font-size: 0.9rem !important;
        font-style: italic;
    }
    
    /* Success/Error message boxes */
    .element-container .stMarkdown .stSuccess {
        background: linear-gradient(135deg, #e8f5e9 0%, #f1f8f4 100%) !important;
        border-left: 4px solid #2e7d32 !important;
        border-radius: 10px !important;
        padding: 1rem !important;
    }
    .element-container .stMarkdown .stError {
        background: linear-gradient(135deg, #ffebee 0%, #fff5f5 100%) !important;
        border-left: 4px solid #c62828 !important;
        border-radius: 10px !important;
        padding: 1rem !important;
    }
    
    /* Links styling */
    a {
        color: #1976d2 !important;
        text-decoration: none !important;
        font-weight: 600 !important;
        transition: all 0.2s ease;
    }
    a:hover {
        color: #0d47a1 !important;
        text-decoration: underline !important;
    }
    
    /* Tooltip improvements */
    [data-testid="stTooltipIcon"] {
        color: #64b5f6 !important;
    }
    
    /* Enhanced hr separator */
    hr {
        border: none !important;
        height: 2px !important;
        background: linear-gradient(90deg, transparent 0%, #bbdefb 50%, transparent 100%) !important;
        margin: 2rem 0 !important;
    }
    </style>
""", unsafe_allow_html=True)


# ---------- CONSTANTS ----------
# Use classification config if available, otherwise fallback to defaults
CONFIG_AVAILABLE = False  # Configuration not available, use fallback
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

# Update color mapping to match the consolidated system
MOTIF_COLORS = MOTIF_CLASS_COLORS  # Use the consolidated colors

PAGES = {
    "Home": "Overview",
    "Upload & Analyze": "Sequence Upload and Motif Analysis",
    "Results": "Analysis Results and Visualization",
    "Download": "Export Data",
    "Documentation": "Scientific Documentation & References"
}

# Remove Entrez setup since it's optional
if BIO_AVAILABLE:
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
EXAMPLE_MULTI_FASTA = """>G4_iMotif_APhilic_Sequence
GGGATTGGGATTGGGATTGGGCCCATCCCTACCCTACCCAAACCCATCCCTACCCTACCCAATTTATTTAAAAA
AAAAAAAAAAAAAAAAAAAAAAGATCGAAAGATCGAAAGATCGAAAGATCGATGCGGCGGCGGCGGCGGCGGCGG
CGGCGGCGAATTCGAATTCGAATTCGAATTCCGCGCGCGCGCGCGCGCGCGAATGCATGCATGCATGCATGCAT
>Z_DNA_RLoop_Complex
CGCGCGCGCGCGCGCGCGCGCGATATATATATATATATATATCGCGCGCGCGCGCGCGCGCGGGGATGGGGATGG
GGATGGGGGGATGGGGATGGGGATGGGTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTGAAA
GAAAAAAGAAAGAAAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAG
>CurvedDNA_SlippedDNA_STR
AAAAAACGTTGCAAAAAACGTTGCAAAAAACGTTGCAAAAAATTTTTTCGAACGTTTTTTCGAACGTTTTTTCGA
ACGCAGCAGCAGCAGCAGCAGCAGCAGCTGCTGCTGCTGCTGCTGCTGCTGATCTGATCTGATCTGATCTGATC
TGATCTGATTCTATTCTATTCTATTCTATTCTATTCTATTCTGGCCCCGGCCCCGGCCCCGGCCCCTGCTGCTG
>Cruciform_Triplex_Mirror
ATGCCCGGGATCGGATCCGATCGAAATTCGATCGGATCCGATCCCGGGCATGAAAGAAAGAAAGAAAGAAAGAAA
GAAAGAAAGAAAAGATCCGGCCGATAGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGTTCCTCCTCCTCC
TCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCGAATTCCGAATTCCGAATTCCGAATTCCGAATTCGAAA
>Multi_iMotif_AC_Sequence  
AAATTTATTTAAATTTAAATTCCCTACCCTACCCTACCCAAAAATCCCTACCCTACCCTACCCGGAATCGATCG
ATCGATCGATCGATCGATCGCCCTACCCTACCCTACCCAAACCCTACCCTACCCTACCCAAAAAAAAAAAAAAAA
AAAAAAAAAAAGATCTAGATCTAGATCTAGATCTAGATCTAGATCTGAAAGAAAGAAAGAAAGAAAGAAAGAAA
"""

# Streamlined session state using consolidated system
for k, v in {
    'seqs': [],
    'names': [],
    'results': [],
    'summary_df': pd.DataFrame(),
    'analysis_status': "Ready",
    'selected_classes': []  # Initialize empty list for motif class selection
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
        <div style='font-family: Inter, system-ui, sans-serif; font-size:1.05rem; color:#263238; 
                    line-height:1.8; padding:2rem; background:white; border-radius:16px; 
                    box-shadow:0 4px 16px rgba(0,0,0,0.08); border: 1px solid rgba(25, 118, 210, 0.1);'>
        <p style='margin-top:0;'><b style='color:#0d47a1; font-size:1.15rem;'>Non-canonical DNA structures</b> play key roles in genome stability, regulation, and evolution.</p>
        <p>This application detects and analyzes <b style='color:#1976d2;'>11 major classes with 22+ subclasses</b> of Non-B DNA motifs in any DNA sequence or multi-FASTA file.</p>
        <p style='margin-bottom:0.8rem;'><b style='color:#0d47a1;'>Motif Classes (11 classes, 22+ subclasses):</b></p>
        <div style='color:#37474f; font-size:0.98rem; line-height:1.9; padding-left:1rem;'>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>1. Curved DNA</b> <span style='color:#546e7a;'>(Global curvature, Local Curvature)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>2. Slipped DNA</b> <span style='color:#546e7a;'>(Direct Repeat, STR)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>3. Cruciform DNA</b> <span style='color:#546e7a;'>(Inverted Repeats)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>4. R-loop</b> <span style='color:#546e7a;'>(R-loop formation sites)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>5. Triplex</b> <span style='color:#546e7a;'>(Triplex, Sticky DNA)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>6. G-Quadruplex Family</b> <span style='color:#546e7a;'>(Multimeric G4, Canonical G4, Relaxed G4, Bulged G4, Bipartite G4, Imperfect G4, G-Triplex intermediate)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>7. i-Motif Family</b> <span style='color:#546e7a;'>(Canonical i-motif, Relaxed i-motif, AC-motif)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>8. Z-DNA</b> <span style='color:#546e7a;'>(Z-DNA, eGZ (Extruded-G) DNA)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>9. A-philic DNA</b> <span style='color:#546e7a;'>(A-philic DNA)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>10. Hybrid</b> <span style='color:#546e7a;'>(dynamic overlaps)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>11. Non-B DNA Clusters</b> <span style='color:#546e7a;'>(dynamic clusters)</span></div>
        </div>
        <p style='margin-top:1.2rem; margin-bottom:0; color:#1976d2; font-weight:600;'>📤 Upload single or multi-FASTA files to begin analysis...</p>
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
                            seqs.append(cur_seq)
                            names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                        cur_name = line.strip().lstrip(">")
                        cur_seq = ""
                    else:
                        cur_seq += line.strip()
                if cur_seq:
                    seqs.append(cur_seq)
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
                            seqs.append(cur_seq)
                            names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                        cur_name = line.strip().lstrip(">")
                        cur_seq = ""
                    else:
                        cur_seq += line.strip()
                if cur_seq:
                    seqs.append(cur_seq)
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
                    parsed_fasta = parse_fasta(EXAMPLE_FASTA)
                    seqs = list(parsed_fasta.values())
                    names = list(parsed_fasta.keys())
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
                                seqs.append(cur_seq)
                                names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                            cur_name = line.strip().lstrip(">")
                            cur_seq = ""
                        else:
                            cur_seq += line.strip()
                    if cur_seq:
                        seqs.append(cur_seq)
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
        
        # Analysis controls simplified
        st.markdown("### ⚙️ Analysis Options")
        
        # Enable all classes by default in consolidated system
        st.info("🔬 **NBDScanner detects all 11 motif classes with 22+ subclasses automatically**")
        
        # Simple options
        col1, col2 = st.columns(2)
        with col1:
            detailed_output = st.checkbox("Detailed Analysis", value=True, 
                                        help="Include comprehensive motif metadata")
        with col2:
            quality_check = st.checkbox("Quality Validation", value=True, 
                                      help="Validate detected motifs")
        
        
        # Hardcoded default overlap handling: always remove overlaps within subclasses
        nonoverlap = True
        overlap_option = "Remove overlaps within subclasses"
        
        
        # ========== RUN ANALYSIS BUTTON ========== 
        if st.button("🔬 Run NBDScanner Analysis", type="primary", use_container_width=True, key="run_motif_analysis_main"):
            # Simplified validation
            if not st.session_state.seqs:
                st.error("❌ Please upload or input sequences before running analysis.")
                st.session_state.analysis_status = "Error"
            else:
                st.session_state.analysis_status = "Running"
                
                # Store analysis parameters in session state for use in download section
                st.session_state.overlap_option_used = overlap_option
                st.session_state.nonoverlap_used = nonoverlap
                
                # Set analysis parameters based on user selections
                # nonoverlap is already set above based on user selection
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
                
                # Enhanced progress tracking with timer
                import time
                
                # Create placeholder for timer and progress
                timer_placeholder = st.empty()
                progress_placeholder = st.empty()
                status_placeholder = st.empty()
                
                start_time = time.time()
                
                try:
                    # Filter which classes to analyze based on selection
                    analysis_classes = st.session_state.selected_classes if st.session_state.selected_classes else None
                    
                    # Run analysis on each sequence
                    all_results = []
                    all_hotspots = []
                    
                    total_bp_processed = 0
                    
                    with progress_placeholder.container():
                        pbar = st.progress(0)
                        
                    for i, (seq, name) in enumerate(zip(st.session_state.seqs, st.session_state.names)):
                        progress = (i + 1) / len(st.session_state.seqs)
                        
                        # Calculate elapsed time
                        elapsed = time.time() - start_time
                        
                        # Update timer display
                        timer_placeholder.markdown(f"""
                        <div style='background: linear-gradient(135deg, #1976d2 0%, #42a5f5 100%); 
                                    border-radius: 12px; padding: 1rem; color: white; text-align: center;
                                    box-shadow: 0 4px 12px rgba(25, 118, 210, 0.3); margin-bottom: 1rem;'>
                            <h3 style='margin: 0; color: white;'>⏱️ Analysis Progress</h3>
                            <h2 style='margin: 0.5rem 0; color: #FFD700; font-size: 2rem;'>{elapsed:.1f}s</h2>
                            <p style='margin: 0; opacity: 0.9;'>Sequence {i+1}/{len(st.session_state.seqs)} | {total_bp_processed:,} bp processed</p>
                        </div>
                        """, unsafe_allow_html=True)
                        
                        # Run the consolidated NBDScanner analysis
                        seq_start = time.time()
                        results = analyze_sequence(seq, name)
                        seq_time = time.time() - seq_start
                        
                        # Ensure all motifs have required fields
                        results = [ensure_subclass(motif) for motif in results]
                        all_results.append(results)
                        
                        total_bp_processed += len(seq)
                        
                        # Calculate processing speed
                        speed = total_bp_processed / elapsed if elapsed > 0 else 0
                        
                        with progress_placeholder.container():
                            pbar.progress(progress, text=f"Analyzed {i+1}/{len(st.session_state.seqs)} sequences")
                        
                        status_placeholder.info(f"✓ {name}: {len(seq):,} bp in {seq_time:.2f}s ({len(seq)/seq_time:.0f} bp/s) - {len(results)} motifs found")
                    
                    # Store results
                    st.session_state.results = all_results
                    
                    # Final timing statistics
                    total_time = time.time() - start_time
                    overall_speed = total_bp_processed / total_time if total_time > 0 else 0
                    
                    # Generate summary
                    summary = []
                    for i, results in enumerate(all_results):
                        seq = st.session_state.seqs[i]
                        stats = get_basic_stats(seq, results)
                        summary.append({
                            'Sequence': st.session_state.names[i],
                            'Length': stats['Length'],
                            'GC Content': f"{stats['GC%']:.1f}%",
                            'Motifs Found': len(results),
                            'Unique Types': len(set(m.get('Type', 'Unknown') for m in results)),
                            'Avg Score': f"{np.mean([m.get('Score', 0) for m in results]):.3f}" if results else "0.000"
                        })
                    
                    st.session_state.summary_df = pd.DataFrame(summary)
                    
                    # Store performance metrics
                    st.session_state.performance_metrics = {
                        'total_time': total_time,
                        'total_bp': total_bp_processed,
                        'speed': overall_speed,
                        'sequences': len(st.session_state.seqs),
                        'total_motifs': sum(len(r) for r in all_results)
                    }
                    
                    # Clear progress displays
                    progress_placeholder.empty()
                    status_placeholder.empty()
                    
                    # Show final success message with performance metrics
                    timer_placeholder.markdown(f"""
                    <div style='background: linear-gradient(135deg, #2e7d32 0%, #4caf50 100%); 
                                border-radius: 12px; padding: 1.5rem; color: white;
                                box-shadow: 0 4px 12px rgba(76, 175, 80, 0.3); margin-bottom: 1rem;'>
                        <h3 style='margin: 0 0 1rem 0; color: white;'>✅ Analysis Complete!</h3>
                        <div style='display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 1rem;'>
                            <div style='text-align: center;'>
                                <h2 style='margin: 0; color: #FFD700; font-size: 1.8rem;'>{total_time:.2f}s</h2>
                                <p style='margin: 0.3rem 0 0 0; opacity: 0.9; font-size: 0.9rem;'>Total Time</p>
                            </div>
                            <div style='text-align: center;'>
                                <h2 style='margin: 0; color: #FFD700; font-size: 1.8rem;'>{total_bp_processed:,}</h2>
                                <p style='margin: 0.3rem 0 0 0; opacity: 0.9; font-size: 0.9rem;'>Base Pairs</p>
                            </div>
                            <div style='text-align: center;'>
                                <h2 style='margin: 0; color: #FFD700; font-size: 1.8rem;'>{overall_speed:,.0f}</h2>
                                <p style='margin: 0.3rem 0 0 0; opacity: 0.9; font-size: 0.9rem;'>bp/second</p>
                            </div>
                            <div style='text-align: center;'>
                                <h2 style='margin: 0; color: #FFD700; font-size: 1.8rem;'>{sum(len(r) for r in all_results)}</h2>
                                <p style='margin: 0.3rem 0 0 0; opacity: 0.9; font-size: 0.9rem;'>Motifs Found</p>
                            </div>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)
                    
                    st.success("Results are available below and in the 'Analysis Results and Visualization' tab.")
                    st.session_state.analysis_status = "Complete"
                    
                except Exception as e:
                    timer_placeholder.empty()
                    progress_placeholder.empty()
                    status_placeholder.empty()
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
        # Performance metrics display if available
        if st.session_state.get('performance_metrics'):
            metrics = st.session_state.performance_metrics
            st.markdown(f"""
            <div style='background: linear-gradient(135deg, #0d47a1 0%, #1976d2 100%); 
                        border-radius: 12px; padding: 1.5rem; color: white; margin-bottom: 2rem;
                        box-shadow: 0 4px 16px rgba(13, 71, 161, 0.25);'>
                <h3 style='margin: 0 0 1rem 0; color: white; text-align: center;'>⚡ Performance Metrics</h3>
                <div style='display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr)); gap: 1rem;'>
                    <div style='text-align: center; background: rgba(255,255,255,0.15); padding: 1rem; border-radius: 8px;'>
                        <h2 style='margin: 0; color: #FFD700; font-size: 1.8rem;'>{metrics['total_time']:.2f}s</h2>
                        <p style='margin: 0.3rem 0 0 0; opacity: 0.9; font-size: 0.85rem;'>Processing Time</p>
                    </div>
                    <div style='text-align: center; background: rgba(255,255,255,0.15); padding: 1rem; border-radius: 8px;'>
                        <h2 style='margin: 0; color: #FFD700; font-size: 1.8rem;'>{metrics['total_bp']:,}</h2>
                        <p style='margin: 0.3rem 0 0 0; opacity: 0.9; font-size: 0.85rem;'>Base Pairs</p>
                    </div>
                    <div style='text-align: center; background: rgba(255,255,255,0.15); padding: 1rem; border-radius: 8px;'>
                        <h2 style='margin: 0; color: #FFD700; font-size: 1.8rem;'>{metrics['speed']:,.0f}</h2>
                        <p style='margin: 0.3rem 0 0 0; opacity: 0.9; font-size: 0.85rem;'>bp/second</p>
                    </div>
                    <div style='text-align: center; background: rgba(255,255,255,0.15); padding: 1rem; border-radius: 8px;'>
                        <h2 style='margin: 0; color: #FFD700; font-size: 1.8rem;'>{metrics['sequences']}</h2>
                        <p style='margin: 0.3rem 0 0 0; opacity: 0.9; font-size: 0.85rem;'>Sequences</p>
                    </div>
                    <div style='text-align: center; background: rgba(255,255,255,0.15); padding: 1rem; border-radius: 8px;'>
                        <h2 style='margin: 0; color: #FFD700; font-size: 1.8rem;'>{metrics['total_motifs']}</h2>
                        <p style='margin: 0.3rem 0 0 0; opacity: 0.9; font-size: 0.85rem;'>Total Motifs</p>
                    </div>
                </div>
            </div>
            """, unsafe_allow_html=True)
        
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
            # Filter motifs for main display (exclude hybrid and cluster)
            filtered_motifs = [m for m in motifs if m.get('Class') not in ['Hybrid', 'Non-B_DNA_Clusters']]
            hybrid_cluster_motifs = [m for m in motifs if m.get('Class') in ['Hybrid', 'Non-B_DNA_Clusters']]
            
            # Create enhanced motifs DataFrame (only non-hybrid/cluster motifs)
            df = pd.DataFrame(filtered_motifs) if filtered_motifs else pd.DataFrame()
            
            # Calculate and display enhanced coverage statistics (using only filtered motifs)
            stats = get_basic_stats(st.session_state.seqs[seq_idx], filtered_motifs)
            
            motif_count = len(filtered_motifs)
            hybrid_cluster_count = len(hybrid_cluster_motifs)
            coverage_pct = stats.get("Motif Coverage %", 0)
            non_b_density = (motif_count / sequence_length * 1000) if sequence_length > 0 else 0
            
            # Enhanced summary card with modern research-quality styling
            st.markdown(f"""
            <div style='background: linear-gradient(135deg, #1976d2 0%, #1565c0 100%); 
                        border-radius: 16px; padding: 2rem; margin: 1.5rem 0; color: white;
                        box-shadow: 0 8px 24px rgba(25, 118, 210, 0.25);'>
                <h3 style='margin: 0 0 1.5rem 0; color: white; text-align: center; font-size: 1.5rem; font-weight: 700;'>
                    🧬 NBDScanner Analysis Results
                </h3>
                <div style='display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); 
                            gap: 1.5rem; margin-top: 1rem;'>
                    <div style='text-align: center; background: rgba(255,255,255,0.15); 
                                padding: 1.2rem; border-radius: 12px; backdrop-filter: blur(10px);'>
                        <h2 style='margin: 0 0 0.5rem 0; color: #FFD700; font-size: 2.2rem; font-weight: 800;'>
                            {stats.get("Coverage%", 0):.2f}%
                        </h2>
                        <p style='margin: 0; font-size: 0.95rem; opacity: 0.95; font-weight: 500;'>
                            Sequence Coverage
                        </p>
                    </div>
                    <div style='text-align: center; background: rgba(255,255,255,0.15); 
                                padding: 1.2rem; border-radius: 12px; backdrop-filter: blur(10px);'>
                        <h2 style='margin: 0 0 0.5rem 0; color: #FFD700; font-size: 2.2rem; font-weight: 800;'>
                            {stats.get("Density", 0):.2f}
                        </h2>
                        <p style='margin: 0; font-size: 0.95rem; opacity: 0.95; font-weight: 500;'>
                            Motif Density<br>(motifs/kb)
                        </p>
                    </div>
                    <div style='text-align: center; background: rgba(255,255,255,0.15); 
                                padding: 1.2rem; border-radius: 12px; backdrop-filter: blur(10px);'>
                        <h2 style='margin: 0 0 0.5rem 0; color: #FFD700; font-size: 2.2rem; font-weight: 800;'>
                            {motif_count}
                        </h2>
                        <p style='margin: 0; font-size: 0.95rem; opacity: 0.95; font-weight: 500;'>
                            Total Motifs
                        </p>
                    </div>
                    <div style='text-align: center; background: rgba(255,255,255,0.15); 
                                padding: 1.2rem; border-radius: 12px; backdrop-filter: blur(10px);'>
                        <h2 style='margin: 0 0 0.5rem 0; color: #FFD700; font-size: 2.2rem; font-weight: 800;'>
                            {sequence_length:,}
                        </h2>
                        <p style='margin: 0; font-size: 0.95rem; opacity: 0.95; font-weight: 500;'>
                            Sequence Length (bp)
                        </p>
                    </div>
                </div>
            </div>
            """, unsafe_allow_html=True)
            
            # Add info about hybrid/cluster motifs being shown separately
            if hybrid_cluster_count > 0:
                st.info(f"ℹ️ {hybrid_cluster_count} Hybrid/Cluster motifs detected. View them in the 'Cluster/Hybrid' tab below.")
            
            # Motif class distribution summary (no score visualization as per requirements)
            if filtered_motifs:
                st.markdown("### 📊 Motif Class Distribution")
                
                # Motif class distribution - show all classes even if 0 (excluding Hybrid and Cluster)
                motif_classes = [m.get('Class', 'Unknown') for m in filtered_motifs]
                
                # Initialize all classes with 0 counts (excluding Hybrid and Cluster)
                filtered_order = [cls for cls in MOTIF_ORDER if cls not in ['Hybrid', 'Non-B DNA Clusters']]
                all_class_counts = {cls: 0 for cls in filtered_order}
                
                # Update with actual counts
                actual_counts = Counter(motif_classes)
                all_class_counts.update(actual_counts)
                
                # Remove 'Unknown' if it exists and has 0 count
                if 'Unknown' in all_class_counts and all_class_counts['Unknown'] == 0:
                    del all_class_counts['Unknown']
                
                fig, ax = plt.subplots(figsize=(10, 6))
                
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
                
                # Create bar chart for class distribution
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
            
            # Column selection for display (remove only Normalized_Score column)
            available_columns = df.columns.tolist()
            # Filter out specific Normalized_Score column
            available_columns = [col for col in available_columns if col != 'Normalized_Score']
            display_columns = st.multiselect(
                "Select columns to display:",
                available_columns,
                default=[col for col in ['Class', 'Subclass', 'Start', 'End', 'Length', 'Score', 'GC Content'] if col in available_columns]
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
            
            # CONSOLIDATED VISUALIZATION SUITE
            st.markdown('<h3>📊 NBDScanner Visualizations</h3>', unsafe_allow_html=True)
            
            # Create tabs for different visualization categories including Cluster/Hybrid tab
            viz_tabs = st.tabs(["📈 Distribution", "🗺️ Coverage Map", "📊 Statistics", "🎯 Interactive", "🔗 Cluster/Hybrid"])
            
            with viz_tabs[0]:  # Distribution
                st.subheader("Motif Distribution Analysis")
                try:
                    fig1 = plot_motif_distribution(filtered_motifs, by='Class', title=f"Motif Classes - {sequence_name}")
                    st.pyplot(fig1)
                    plt.close(fig1)
                    
                    fig2 = plot_motif_distribution(filtered_motifs, by='Subclass', title=f"Motif Subclasses - {sequence_name}")
                    st.pyplot(fig2) 
                    plt.close(fig2)
                    
                    # Pie chart
                    fig3 = plot_nested_pie_chart(filtered_motifs, title=f"Class-Subclass Distribution - {sequence_name}")
                    st.pyplot(fig3)
                    plt.close(fig3)
                except Exception as e:
                    st.error(f"Error generating distribution plots: {e}")
            
            with viz_tabs[1]:  # Coverage Map
                st.subheader("Sequence Coverage Analysis")
                try:
                    fig4 = plot_coverage_map(filtered_motifs, sequence_length, title=f"Motif Coverage - {sequence_name}")
                    st.pyplot(fig4)
                    plt.close(fig4)
                except Exception as e:
                    st.error(f"Error generating coverage map: {e}")
            
            with viz_tabs[2]:  # Statistics  
                st.subheader("Statistical Analysis")
                try:
                    # Only show length distribution, not score distribution as per requirements
                    fig6 = plot_length_distribution(filtered_motifs, by_class=True, title="Length Distribution by Class") 
                    st.pyplot(fig6)
                    plt.close(fig6)
                except Exception as e:
                    st.error(f"Error generating statistical plots: {e}")
            
            with viz_tabs[3]:  # Interactive
                st.subheader("Interactive Analysis")
                try:
                    from visualization import create_interactive_coverage_plot
                    interactive_fig = create_interactive_coverage_plot(filtered_motifs, sequence_length, 
                                                                     title=f"Interactive Motif Browser - {sequence_name}")
                    if hasattr(interactive_fig, 'show'):  # Plotly figure
                        st.plotly_chart(interactive_fig, use_container_width=True)
                    else:  # Matplotlib fallback
                        st.pyplot(interactive_fig)
                        plt.close(interactive_fig)
                except Exception as e:
                    st.error(f"Error generating interactive plots: {e}")
                    st.info("Interactive plots require plotly. Install with: pip install plotly")
                    
                    # Simple fallback visualization for error cases
                    st.markdown('<h4>📊 Fallback Visualization</h4>', unsafe_allow_html=True)
                    st.info("Using basic fallback visualizations since comprehensive analysis failed.")
                    
                    # Basic motif count chart
                    st.markdown('<span style="font-family:Montserrat,Arial;font-size:1.11rem;"><b>Motif Class Distribution</b></span>', unsafe_allow_html=True)
                    if filtered_motifs:
                        df_viz = pd.DataFrame(filtered_motifs)
                        fig, ax = plt.subplots(figsize=(10,6))
                        class_counts = df_viz['Class'].value_counts()
                        ax.barh(class_counts.index, class_counts.values)
                        ax.set_xlabel("Motif Count")
                        ax.set_title("Basic Motif Class Distribution")
                        st.pyplot(fig)
                        plt.close(fig)
                        
                        # Basic motif map
                        st.markdown('<span style="font-family:Montserrat,Arial;font-size:1.11rem;"><b>Motif Position Map</b></span>', unsafe_allow_html=True)
                        fig, ax = plt.subplots(figsize=(12,4))
                        for i, (_, row) in enumerate(df_viz.iterrows()):
                            ax.plot([row['Start'], row['End']], [i, i], lw=3, alpha=0.7)
                        ax.set_xlabel("Sequence Position (bp)")
                        ax.set_ylabel("Motif Index")
                        ax.set_title("Basic Motif Position Map")
                        st.pyplot(fig)
                        plt.close(fig)
            
            with viz_tabs[4]:  # Cluster/Hybrid
                st.subheader("🔗 Hybrid and Cluster Motifs")
                
                if not hybrid_cluster_motifs:
                    st.info("No Hybrid or Cluster motifs detected for this sequence.")
                else:
                    st.markdown(f"""
                    <div style='background: #f0f8ff; border-left: 4px solid #1e88e5; padding: 15px; border-radius: 5px; margin: 10px 0;'>
                    <h4 style='margin-top: 0; color: #1565c0;'>About Hybrid and Cluster Motifs</h4>
                    <p><b>Hybrid Motifs:</b> Regions where different non-B DNA classes overlap (30-70% overlap between classes).</p>
                    <p><b>Cluster Motifs:</b> High-density regions containing multiple non-B DNA motifs from different classes.</p>
                    <p>These motifs represent complex genomic regions with potential biological significance.</p>
                    </div>
                    """, unsafe_allow_html=True)
                    
                    # Create DataFrame for hybrid/cluster motifs
                    hc_df = pd.DataFrame(hybrid_cluster_motifs)
                    
                    # Summary statistics
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        hybrid_count = len([m for m in hybrid_cluster_motifs if m.get('Class') == 'Hybrid'])
                        st.metric("Hybrid Motifs", hybrid_count)
                    with col2:
                        cluster_count = len([m for m in hybrid_cluster_motifs if m.get('Class') == 'Non-B_DNA_Clusters'])
                        st.metric("Cluster Motifs", cluster_count)
                    with col3:
                        avg_length = int(hc_df['Length'].mean()) if 'Length' in hc_df.columns else 0
                        st.metric("Avg Length (bp)", avg_length)
                    
                    # Detailed table
                    st.markdown("### 📋 Detailed Cluster/Hybrid Table")
                    display_cols = ['Class', 'Subclass', 'Start', 'End', 'Length', 'Score']
                    if 'Component_Classes' in hc_df.columns:
                        display_cols.append('Component_Classes')
                    if 'Motif_Count' in hc_df.columns:
                        display_cols.append('Motif_Count')
                    if 'Class_Diversity' in hc_df.columns:
                        display_cols.append('Class_Diversity')
                    
                    # Filter to only show available columns
                    available_display_cols = [col for col in display_cols if col in hc_df.columns]
                    display_hc_df = hc_df[available_display_cols].copy()
                    display_hc_df.columns = [col.replace('_', ' ') for col in display_hc_df.columns]
                    st.dataframe(display_hc_df, use_container_width=True, height=300)
                    
                    # Visualizations for hybrid/cluster
                    st.markdown("### 📊 Visualizations")
                    
                    viz_col1, viz_col2 = st.columns(2)
                    
                    with viz_col1:
                        # Position map
                        st.markdown("**Position Map**")
                        fig, ax = plt.subplots(figsize=(10, 4))
                        colors_map = {'Hybrid': '#ff6b6b', 'Non-B_DNA_Clusters': '#4ecdc4'}
                        for i, motif in enumerate(hybrid_cluster_motifs):
                            color = colors_map.get(motif.get('Class'), '#95a5a6')
                            ax.barh(i, motif['End'] - motif['Start'], left=motif['Start'], 
                                   height=0.8, color=color, alpha=0.7, label=motif.get('Class'))
                        ax.set_xlabel("Sequence Position (bp)")
                        ax.set_ylabel("Motif Index")
                        ax.set_title("Hybrid/Cluster Position Map")
                        # Remove duplicate labels
                        handles, labels = ax.get_legend_handles_labels()
                        by_label = dict(zip(labels, handles))
                        ax.legend(by_label.values(), by_label.keys())
                        plt.tight_layout()
                        st.pyplot(fig)
                        plt.close(fig)
                    
                    with viz_col2:
                        # Raw score distribution (no normalized scores as per requirements)
                        st.markdown("**Raw Score Distribution**")
                        fig, ax = plt.subplots(figsize=(10, 4))
                        hybrid_scores = [m.get('Score', 0) for m in hybrid_cluster_motifs if m.get('Class') == 'Hybrid']
                        cluster_scores = [m.get('Score', 0) for m in hybrid_cluster_motifs if m.get('Class') == 'Non-B_DNA_Clusters']
                        
                        if hybrid_scores:
                            ax.hist(hybrid_scores, bins=10, alpha=0.6, color='#ff6b6b', label='Hybrid')
                        if cluster_scores:
                            ax.hist(cluster_scores, bins=10, alpha=0.6, color='#4ecdc4', label='Cluster')
                        
                        ax.set_xlabel("Raw Score")
                        ax.set_ylabel("Frequency")
                        ax.set_title("Raw Score Distribution")
                        ax.legend()
                        plt.tight_layout()
                        st.pyplot(fig)
                        plt.close(fig)

# ---------- DOWNLOAD ----------
with tab_pages["Download"]:
    st.header("Export Data")
    if not st.session_state.results:
        st.info("No results available to download.")
    else:
        st.markdown("### 📊 Export Options")
        
        # Show analysis configuration used
        overlap_option_used = st.session_state.get('overlap_option_used', 'Remove overlaps within subclasses (default)')
        st.info(f"""
        **Analysis Configuration Used:**
        • Overlap Handling: {overlap_option_used}
        • Motif Classes: {len(st.session_state.get('selected_classes', []))} classes selected
        • Total Motifs Found: {sum(len(motifs) for motifs in st.session_state.results)}
        """)
        
        # Export configuration
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**📊 Export Configuration**")
            
        with col2:
            include_sequences = st.checkbox("Include Full Sequences", value=True,
                                           help="Include full motif sequences in export")
        
        # Prepare data for export using consolidated utilities (exclude Hybrid and Cluster)
        all_motifs = []
        hybrid_cluster_motifs_export = []
        for i, motifs in enumerate(st.session_state.results):
            for m in motifs:
                export_motif = m.copy()
                if 'Sequence_Name' not in export_motif:
                    export_motif['Sequence_Name'] = st.session_state.names[i]
                
                # Separate hybrid/cluster from regular motifs
                if export_motif.get('Class') in ['Hybrid', 'Non-B_DNA_Clusters']:
                    hybrid_cluster_motifs_export.append(export_motif)
                else:
                    all_motifs.append(export_motif)
        
        # Info about excluded motifs
        if hybrid_cluster_motifs_export:
            st.info(f"ℹ️ {len(hybrid_cluster_motifs_export)} Hybrid/Cluster motifs are excluded from downloads. These are shown only in the Cluster/Hybrid visualization tab.")
        
        # Display preview
        if all_motifs:
            df_preview = export_results_to_dataframe(all_motifs).head(10)
            st.markdown("### 👀 Export Preview")
            st.dataframe(df_preview, use_container_width=True)
            st.caption(f"Showing first 10 of {len(all_motifs)} total records (Hybrid/Cluster motifs excluded)")
        
        # Export buttons using consolidated functions
        st.markdown("### 💾 Download Files")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            # CSV Export
            if all_motifs:
                csv_data = export_to_csv(all_motifs)
                st.download_button(
                    "📄 Download CSV", 
                    data=csv_data.encode('utf-8'), 
                    file_name="nbdscanner_results.csv", 
                    mime="text/csv",
                    use_container_width=True,
                    help="Comma-separated values format"
                )
        
        with col2:
            # JSON Export  
            if all_motifs:
                json_data = export_to_json(all_motifs, pretty=True)
                st.download_button(
                    "📊 Download JSON", 
                    data=json_data.encode('utf-8'), 
                    file_name="nbdscanner_results.json", 
                    mime="application/json",
                    use_container_width=True,
                    help="JSON format with metadata"
                )
        
        with col3:
            # BED Export
            if all_motifs and st.session_state.names:
                bed_data = export_to_bed(all_motifs, st.session_state.names[0])
                st.download_button(
                    "🧬 Download BED", 
                    data=bed_data.encode('utf-8'), 
                    file_name="nbdscanner_results.bed", 
                    mime="text/plain",
                    use_container_width=True,
                    help="BED format for genome browsers"
                )
            
            # Config Export - always available when results exist
            if all_motifs:
                import json
                # Create a basic configuration summary
                config_summary = {
                    "analysis_info": {
                        "tool": "NBDScanner",
                        "version": "2024.1",
                        "motif_classes_detected": list(set(m.get('Class', 'Unknown') for m in all_motifs)),
                        "total_sequences_analyzed": len(st.session_state.names) if hasattr(st.session_state, 'names') else 0,
                        "total_motifs_found": len(all_motifs)
                    },
                    "detection_parameters": {
                        "all_classes_enabled": True,
                        "consolidated_system": True
                    }
                }
                config_json = json.dumps(config_summary, indent=2)
                st.download_button(
                    "⚙️ Download Config", 
                    data=config_json, 
                    file_name="analysis_configuration.json", 
                    mime="application/json",
                    use_container_width=True,
                    help="Analysis configuration and metadata"
                )
        
        # Prepare data for genome browser exports
        all_seq_data = []
        for i, motifs in enumerate(st.session_state.results):
            all_seq_data.append({
                'sequence_name': st.session_state.names[i],
                'motifs': motifs
            })


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
