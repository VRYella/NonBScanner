"""
Non-B DNA Scanner - Streamlit Application Interface
===================================================

| Interface Section | Function                                          |
|-------------------|---------------------------------------------------|
| Home Page         | Welcome, overview, and tool introduction         |
| Input Methods     | Text/FASTA input, file upload, NCBI fetch        |
| Analysis Config   | Motif class selection and parameter tuning       |
| Results Display   | Tabular results, statistics, visualizations      |
| Export Options    | CSV, Excel, BED format downloads                 |

Main GUI application for Non-B DNA structure prediction and analysis
"""

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import io
import numpy as np
from collections import Counter

from utils import parse_fasta, gc_content, reverse_complement, wrap, get_basic_stats, validate_dna_sequence
from detection_scoring import detect_motifs_in_sequence, analyze_multiple_sequences, calculate_motif_statistics
from registires import CLASS_DEFINITIONS, DEFAULT_SELECTED_CLASSES, get_motif_limits

# =============================================================================
# PAGE CONFIGURATION
# =============================================================================
st.set_page_config(
    page_title="Non-B DNA Scanner",
    layout="wide",
    page_icon="🧬",
    menu_items={'About': "Non-B DNA Scanner | Developed by Dr. Venkata Rajesh Yella"}
)

# =============================================================================
# SESSION STATE INITIALIZATION  
# =============================================================================
if 'sequence' not in st.session_state:
    st.session_state.sequence = ""
if 'sequence_name' not in st.session_state:
    st.session_state.sequence_name = "Input_Sequence"
if 'motifs' not in st.session_state:
    st.session_state.motifs = []
if 'analysis_complete' not in st.session_state:
    st.session_state.analysis_complete = False
if 'selected_classes' not in st.session_state:
    st.session_state.selected_classes = DEFAULT_SELECTED_CLASSES.copy()

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def ensure_subclass(motif):
    """Ensure every motif has a valid Subclass field"""
    if isinstance(motif, dict):
        if 'Subclass' not in motif or motif['Subclass'] is None:
            motif['Subclass'] = motif.get('Subtype', 'Other')
        return motif
    else:
        return {'Subclass': 'Other', 'Motif': motif}

def format_motifs_dataframe(motifs):
    """Format motifs for display in DataFrame with proper column names"""
    if not motifs:
        return pd.DataFrame()
    
    # Ensure all motifs have required fields
    formatted_motifs = [ensure_subclass(motif) for motif in motifs]
    df = pd.DataFrame(formatted_motifs)
    
    # Column mapping for better display
    column_mapping = {
        'Class': 'Motif Class',
        'Subclass': 'Subclass',
        'Start': 'Start Position', 
        'End': 'End Position',
        'Length': 'Length (bp)',
        'Normalized_Score': 'Confidence Score',
        'Raw_Score': 'Raw Score',
        'Scoring_Method': 'Scoring Method',
        'Sequence_Name': 'Sequence Name'
    }
    
    # Rename columns and select relevant ones
    display_columns = [col for col in column_mapping.keys() if col in df.columns]
    df_display = df[display_columns].copy()
    df_display.rename(columns=column_mapping, inplace=True)
    
    # Round numeric columns
    numeric_columns = ['Confidence Score', 'Raw Score']
    for col in numeric_columns:
        if col in df_display.columns:
            df_display[col] = df_display[col].round(3)
    
    return df_display

def get_biological_significance(class_id: int) -> str:
    """Get biological significance description for each structure class"""
    significance_map = {
        1: "Gene regulation, nucleosome positioning",
        2: "Replication slippage, repeat expansion diseases", 
        3: "Recombination hotspots, chromosomal translocations",
        4: "Transcription regulation, genome instability",
        5: "Antigene therapy, DNA damage",
        6: "Telomere maintenance, oncogene regulation",
        7: "pH sensing, gene expression control",
        8: "Transcription regulation, immune recognition",
        9: "Protein-DNA interactions, A-tract formation, transcription regulation"
    }
    return significance_map.get(class_id, "Unknown biological role")

# =============================================================================
# MAIN APPLICATION INTERFACE
# =============================================================================

def main():
    """Main application function"""
    
    # Application header with logo
    col1, col2 = st.columns([1, 4])
    with col1:
        try:
            st.image("nbdcircle.png", width=120)
        except:
            st.markdown("🧬", unsafe_allow_html=True)
    
    with col2:
        st.title("Non-B DNA Scanner")
        st.markdown("**Advanced Detection of Non-Canonical DNA Structures**")
        
    # Navigation tabs
    tabs = st.tabs(["🏠 Home", "📝 Input & Analysis", "📊 Results", "💾 Export"])
    
    # =============================================================================
    # HOME TAB
    # =============================================================================
    with tabs[0]:
        st.markdown("""
        ## Welcome to Non-B DNA Scanner
        
        This tool identifies and analyzes non-canonical DNA structures that deviate from 
        the standard B-form double helix. These structures play important roles in gene 
        regulation, DNA replication, and genome stability.
        
        ### Detected Structure Types
        """)
        
        # Structure types table
        structure_data = []
        for class_id, info in CLASS_DEFINITIONS.items():
            structure_data.append({
                'Class': class_id,
                'Structure Type': info['name'],
                'Description': info['description'],
                'Biological Significance': get_biological_significance(class_id)
            })
        
        df_structures = pd.DataFrame(structure_data)
        st.dataframe(df_structures, use_container_width=True)
        
        st.markdown("""
        ### Key Features
        
        | Feature | Description |
        |---------|-------------|
        | 9 Structure Types | G-quadruplex, i-motif, Z-DNA, curved DNA, triplex, cruciform, R-loops, slipped DNA, A-philic DNA |
        | Scientific Scoring | Literature-based algorithms (G4Hunter, Z-seeker, etc.) |
        | Multiple Input Methods | Text input, file upload |
        | Advanced Filtering | Score thresholds, overlap removal, length constraints |
        | Export Options | CSV, Excel, BED formats for further analysis |
        """)
    
    # =============================================================================
    # INPUT & ANALYSIS TAB
    # =============================================================================
    with tabs[1]:
        st.header("Sequence Input & Analysis Configuration")
        
        # Input method selection
        input_method = st.radio(
            "Choose input method:",
            ["Direct Text Input", "File Upload"]
        )
        
        sequence = ""
        sequence_name = "Input_Sequence"
        
        # Direct text input
        if input_method == "Direct Text Input":
            sequence_input = st.text_area(
                "Enter DNA sequence (FASTA format or plain sequence):",
                height=150,
                placeholder="Paste your DNA sequence here..."
            )
            
            if sequence_input:
                sequence = parse_fasta(sequence_input)
                sequence_name = st.text_input("Sequence name:", value="User_Sequence")
        
        # File upload
        elif input_method == "File Upload":
            uploaded_file = st.file_uploader(
                "Choose a FASTA file",
                type=['fasta', 'fa', 'txt'],
                help="Upload a FASTA format file containing DNA sequences"
            )
            
            if uploaded_file is not None:
                content = uploaded_file.read().decode('utf-8')
                sequence = parse_fasta(content)
                sequence_name = uploaded_file.name.split('.')[0]
        
        # Sequence validation and display
        if sequence:
            is_valid, error_msg = validate_dna_sequence(sequence)
            
            if is_valid:
                st.success(f"✅ Valid DNA sequence: {len(sequence):,} bp")
                
                # Sequence preview
                with st.expander("Sequence Preview"):
                    st.text(wrap(sequence[:500] + ("..." if len(sequence) > 500 else "")))
                
                # Basic statistics
                basic_stats = get_basic_stats(sequence)
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.metric("Length", f"{basic_stats['Length']:,} bp")
                with col2:
                    st.metric("GC Content", f"{basic_stats['GC_Content']:.1f}%")
                with col3:
                    st.metric("AT Content", f"{basic_stats['AT_Content']:.1f}%")
                
                # Analysis configuration
                st.subheader("Analysis Configuration")
                
                # Motif class selection
                selected_classes = st.multiselect(
                    "Select motif classes to analyze:",
                    options=list(CLASS_DEFINITIONS.keys()),
                    default=DEFAULT_SELECTED_CLASSES,
                    format_func=lambda x: f"Class {x}: {CLASS_DEFINITIONS[x]['name']}"
                )
                
                # Advanced parameters
                with st.expander("Advanced Parameters"):
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        min_score = st.slider("Minimum confidence score", 0.0, 1.0, 0.3, 0.1)
                        overlap_threshold = st.slider("Overlap removal threshold", 0.0, 1.0, 0.5, 0.1)
                    
                    with col2:
                        max_gap = st.number_input("Maximum gap for merging (bp)", 1, 100, 10)
                        remove_overlaps = st.checkbox("Remove overlapping motifs", True)
                
                # Analysis button
                if st.button("🔍 Run Analysis", type="primary"):
                    if selected_classes:
                        with st.spinner("Analyzing sequence for Non-B DNA structures..."):
                            # Configure analysis filters
                            filters = {
                                'remove_overlaps': remove_overlaps,
                                'overlap_threshold': overlap_threshold,
                                'score_filter': True,
                                'min_score': min_score,
                                'merge_nearby': True,
                                'max_gap': max_gap
                            }
                            
                            # Run analysis
                            sequences = {sequence_name: sequence}
                            results = analyze_multiple_sequences(sequences, filters=filters)
                            motifs = results.get(sequence_name, [])
                            
                            # Filter by selected classes
                            motifs = [m for m in motifs if m.get('Class') in selected_classes]
                            
                            # Store results in session state
                            st.session_state.sequence = sequence
                            st.session_state.sequence_name = sequence_name
                            st.session_state.motifs = motifs
                            st.session_state.analysis_complete = True
                            
                            st.success(f"Analysis complete! Found {len(motifs)} Non-B DNA motifs.")
                    else:
                        st.warning("Please select at least one motif class for analysis.")
            else:
                st.error(f"❌ Invalid sequence: {error_msg}")
    
    # =============================================================================
    # RESULTS TAB
    # =============================================================================
    with tabs[2]:
        st.header("Analysis Results")
        
        if st.session_state.analysis_complete and st.session_state.motifs:
            motifs = st.session_state.motifs
            
            # Summary statistics
            st.subheader("Summary Statistics")
            stats = calculate_motif_statistics(motifs)
            
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Total Motifs", stats['Total_Motifs'])
            with col2:
                st.metric("Classes Found", stats['Classes_Found'])
            with col3:
                st.metric("Average Score", f"{stats['Average_Score']:.3f}")
            with col4:
                st.metric("Average Length", f"{stats['Average_Length']:.1f} bp")
            
            # Class distribution
            st.subheader("Class Distribution")
            class_dist = stats['Class_Distribution']
            df_dist = pd.DataFrame(list(class_dist.items()), columns=['Motif Class', 'Count'])
            
            col1, col2 = st.columns([2, 1])
            with col1:
                st.bar_chart(df_dist.set_index('Motif Class'))
            with col2:
                st.dataframe(df_dist, use_container_width=True)
            
            # Detailed results table
            st.subheader("Detailed Results")
            df_motifs = format_motifs_dataframe(motifs)
            
            if not df_motifs.empty:
                # Filter options
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    class_filter = st.multiselect(
                        "Filter by class:",
                        options=df_motifs['Motif Class'].unique() if 'Motif Class' in df_motifs.columns else [],
                        default=df_motifs['Motif Class'].unique() if 'Motif Class' in df_motifs.columns else []
                    )
                
                with col2:
                    if 'Confidence Score' in df_motifs.columns:
                        min_score_filter = st.slider(
                            "Minimum score:",
                            float(df_motifs['Confidence Score'].min()),
                            float(df_motifs['Confidence Score'].max()),
                            float(df_motifs['Confidence Score'].min())
                        )
                    else:
                        min_score_filter = 0.0
                
                with col3:
                    if 'Length (bp)' in df_motifs.columns:
                        min_length = st.number_input(
                            "Minimum length (bp):",
                            int(df_motifs['Length (bp)'].min()),
                            int(df_motifs['Length (bp)'].max()),
                            int(df_motifs['Length (bp)'].min())
                        )
                    else:
                        min_length = 0
                
                # Apply filters
                filtered_df = df_motifs.copy()
                if class_filter and 'Motif Class' in filtered_df.columns:
                    filtered_df = filtered_df[filtered_df['Motif Class'].isin(class_filter)]
                if 'Confidence Score' in filtered_df.columns:
                    filtered_df = filtered_df[filtered_df['Confidence Score'] >= min_score_filter]
                if 'Length (bp)' in filtered_df.columns:
                    filtered_df = filtered_df[filtered_df['Length (bp)'] >= min_length]
                
                # Display filtered results
                st.dataframe(filtered_df, use_container_width=True)
                st.caption(f"Showing {len(filtered_df)} of {len(df_motifs)} total motifs")
            else:
                st.info("No motifs to display")
        
        else:
            st.info("No analysis results available. Please run an analysis first.")
    
    # =============================================================================
    # EXPORT TAB
    # =============================================================================
    with tabs[3]:
        st.header("Export Results")
        
        if st.session_state.analysis_complete and st.session_state.motifs:
            motifs = st.session_state.motifs
            df_motifs = format_motifs_dataframe(motifs)
            
            st.subheader("Export Options")
            
            # Preview data
            with st.expander("Data Preview"):
                st.dataframe(df_motifs.head(10), use_container_width=True)
                st.caption(f"Showing first 10 of {len(df_motifs)} total records")
            
            # Export buttons
            col1, col2, col3 = st.columns(3)
            
            with col1:
                # CSV export
                csv_data = df_motifs.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="📄 Download CSV",
                    data=csv_data,
                    file_name=f"{st.session_state.sequence_name}_motifs.csv",
                    mime="text/csv",
                    use_container_width=True
                )
            
            with col2:
                # Excel export
                excel_buffer = io.BytesIO()
                with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
                    df_motifs.to_excel(writer, sheet_name='Motifs', index=False)
                    
                    # Add summary sheet
                    stats = calculate_motif_statistics(motifs)
                    df_summary = pd.DataFrame([stats])
                    df_summary.to_excel(writer, sheet_name='Summary', index=False)
                
                st.download_button(
                    label="📊 Download Excel",
                    data=excel_buffer.getvalue(),
                    file_name=f"{st.session_state.sequence_name}_motifs.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    use_container_width=True
                )
            
            with col3:
                # BED format export
                bed_lines = ["track name=\"NonB_DNA_Motifs\" description=\"Non-B DNA Structures\"\n"]
                for motif in motifs:
                    chrom = st.session_state.sequence_name
                    start = motif.get('Start', 1) - 1  # Convert to 0-based
                    end = motif.get('End', 1)
                    name = f"Class{motif.get('Class', 0)}_{motif.get('Subclass', 'Unknown')}"
                    score = min(1000, int(motif.get('Normalized_Score', 0) * 1000))
                    strand = '+'
                    
                    bed_lines.append(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
                
                bed_data = ''.join(bed_lines).encode('utf-8')
                st.download_button(
                    label="🧬 Download BED",
                    data=bed_data,
                    file_name=f"{st.session_state.sequence_name}_motifs.bed",
                    mime="text/plain",
                    use_container_width=True
                )
            
            # Export statistics
            st.subheader("Analysis Summary Report")
            
            summary_report = f"""
# Non-B DNA Analysis Report

**Sequence:** {st.session_state.sequence_name}
**Length:** {len(st.session_state.sequence):,} bp
**Analysis Date:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

## Summary Statistics
- **Total Motifs Found:** {len(motifs)}
- **Classes Detected:** {len(set(m.get('Class', 0) for m in motifs))}
- **Average Confidence Score:** {np.mean([m.get('Normalized_Score', 0) for m in motifs]):.3f}
- **Average Length:** {np.mean([m.get('Length', 0) for m in motifs]):.1f} bp

## Class Distribution
"""
            
            class_counts = Counter(m.get('Class', 0) for m in motifs)
            for class_id, count in sorted(class_counts.items()):
                class_name = CLASS_DEFINITIONS.get(class_id, {}).get('name', f'Class {class_id}')
                summary_report += f"- **{class_name}:** {count} motifs\n"
            
            st.markdown(summary_report)
            
            # Download report
            st.download_button(
                label="📋 Download Analysis Report",
                data=summary_report.encode('utf-8'),
                file_name=f"{st.session_state.sequence_name}_analysis_report.md",
                mime="text/markdown"
            )
        
        else:
            st.info("No results available for export. Please run an analysis first.")

if __name__ == "__main__":
    main()