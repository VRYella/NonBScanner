"""
NBDFinder I/O Module
==================

Input/Output functionality for NBDFinder including:
- FASTA file handling with memory mapping
- Export utilities for various formats (BED, CSV, Parquet, GFF)
- Pydantic schemas for type-safe data models
"""

from .fasta import (
    FastaReader, parse_fasta_string, write_fasta, extract_sequence_regions
)

from .writers import (
    export_to_bed, export_to_csv, export_to_gff3, export_to_parquet,
    export_to_json, create_track_hub
)

from .schemas import (
    MotifClass, Strand, MotifRecord, SequenceInfo, AnalysisConfig,
    AnalysisResults, ExportConfig, BatchAnalysisRequest, PerformanceMetrics,
    create_motif_record, create_sequence_info
)

__all__ = [
    # FASTA handling
    'FastaReader',
    'parse_fasta_string', 
    'write_fasta',
    'extract_sequence_regions',
    
    # Export functions
    'export_to_bed',
    'export_to_csv',
    'export_to_gff3',
    'export_to_parquet',
    'export_to_json',
    'create_track_hub',
    
    # Data models
    'MotifClass',
    'Strand',
    'MotifRecord',
    'SequenceInfo',
    'AnalysisConfig',
    'AnalysisResults',
    'ExportConfig',
    'BatchAnalysisRequest',
    'PerformanceMetrics',
    
    # Factory functions
    'create_motif_record',
    'create_sequence_info'
]