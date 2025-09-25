"""
Pydantic Models for NBDFinder Data Structures
============================================

Type-safe data models for motif records, configuration, and results
using Pydantic for validation and serialization.
"""

from pydantic import BaseModel, Field, validator, ConfigDict
from typing import List, Dict, Any, Optional, Union, Literal
from datetime import datetime
from enum import Enum

class MotifClass(str, Enum):
    """Enumeration of motif classes."""
    CURVED_DNA = "Curved_DNA"
    SLIPPED_DNA = "Slipped_DNA" 
    CRUCIFORM = "Cruciform"
    R_LOOP = "R-Loop"
    TRIPLEX = "Triplex"
    G_QUADRUPLEX = "G-Quadruplex"
    I_MOTIF = "i-Motif"
    Z_DNA = "Z-DNA"
    HYBRID = "Hybrid"
    CLUSTER = "Cluster"

class Strand(str, Enum):
    """DNA strand orientation."""
    PLUS = "+"
    MINUS = "-"
    BOTH = "±"

class MotifRecord(BaseModel):
    """
    Core motif record with standardized fields.
    """
    model_config = ConfigDict(str_strip_whitespace=True)
    
    # Core identification
    motif_class: MotifClass = Field(..., description="Main motif class")
    subclass: str = Field(..., description="Specific motif subclass")
    
    # Genomic coordinates (1-based, inclusive)
    chromosome: Optional[str] = Field(None, description="Chromosome/sequence name")
    start: int = Field(..., ge=1, description="Start position (1-based)")
    end: int = Field(..., ge=1, description="End position (1-based, inclusive)")
    strand: Strand = Field(Strand.PLUS, description="Strand orientation")
    
    # Sequence information
    sequence: str = Field(..., min_length=1, description="Motif sequence")
    length: Optional[int] = Field(None, ge=1, description="Motif length in bp")
    
    # Scoring information
    score: float = Field(0.0, ge=0.0, description="Raw motif score")
    normalized_score: Optional[float] = Field(None, ge=0.0, le=1.0, description="Normalized score (0-1)")
    score_method: str = Field("default", description="Scoring algorithm used")
    
    # Analysis metadata
    detection_method: str = Field("regex", description="Detection method")
    pattern_id: Optional[int] = Field(None, description="Pattern ID from registry")
    confidence: Optional[float] = Field(None, ge=0.0, le=1.0, description="Detection confidence")
    
    # Additional properties
    properties: Dict[str, Any] = Field(default_factory=dict, description="Additional motif properties")
    
    @validator('length', always=True)
    def calculate_length(cls, v, values):
        """Calculate length from coordinates if not provided."""
        if v is None and 'start' in values and 'end' in values:
            return values['end'] - values['start'] + 1
        return v
    
    @validator('sequence')
    def validate_dna_sequence(cls, v):
        """Validate DNA sequence contains only valid bases."""
        valid_bases = set('ATGCN')
        if not all(base in valid_bases for base in v.upper()):
            raise ValueError("Sequence must contain only A, T, G, C, N characters")
        return v.upper()
    
    @validator('end')
    def end_after_start(cls, v, values):
        """Ensure end position is after start position."""
        if 'start' in values and v < values['start']:
            raise ValueError("End position must be >= start position")
        return v

class SequenceInfo(BaseModel):
    """Information about analyzed sequence."""
    model_config = ConfigDict(str_strip_whitespace=True)
    
    name: str = Field(..., description="Sequence name/identifier")
    length: int = Field(..., ge=1, description="Sequence length in bp")
    gc_content: float = Field(..., ge=0.0, le=100.0, description="GC content percentage")
    
    # Optional metadata
    description: Optional[str] = Field(None, description="Sequence description")
    organism: Optional[str] = Field(None, description="Source organism")
    chromosome: Optional[str] = Field(None, description="Chromosome identifier")
    
    # Analysis statistics
    total_motifs: int = Field(0, ge=0, description="Total motifs detected")
    motif_coverage: float = Field(0.0, ge=0.0, le=100.0, description="Percentage covered by motifs")
    
class AnalysisConfig(BaseModel):
    """Configuration for motif analysis."""
    model_config = ConfigDict(str_strip_whitespace=True)
    
    # Detection parameters
    chunk_size: int = Field(100000, ge=1000, description="Chunk size for large sequences")
    overlap_size: int = Field(1000, ge=0, description="Overlap between chunks")
    max_workers: Optional[int] = Field(None, ge=1, description="Maximum parallel workers")
    
    # Filtering thresholds
    min_score_threshold: float = Field(0.1, ge=0.0, description="Minimum score threshold")
    min_motif_length: int = Field(10, ge=1, description="Minimum motif length")
    max_motif_length: int = Field(1000, ge=1, description="Maximum motif length")
    
    # Post-processing options
    remove_overlaps: bool = Field(True, description="Remove overlapping motifs")
    merge_nearby: bool = Field(False, description="Merge nearby motifs")
    merge_distance: int = Field(10, ge=0, description="Maximum distance for merging")
    
    # Algorithm selection
    enabled_classes: List[MotifClass] = Field(
        default_factory=lambda: list(MotifClass),
        description="Motif classes to detect"
    )
    scoring_methods: Dict[str, str] = Field(
        default_factory=dict,
        description="Scoring method per motif class"
    )

class AnalysisResults(BaseModel):
    """Complete analysis results."""
    model_config = ConfigDict(str_strip_whitespace=True)
    
    # Metadata
    analysis_id: str = Field(..., description="Unique analysis identifier")
    timestamp: datetime = Field(default_factory=datetime.now, description="Analysis timestamp")
    version: str = Field("2.0.0", description="NBDFinder version")
    
    # Input information
    sequence_info: SequenceInfo = Field(..., description="Analyzed sequence information")
    config: AnalysisConfig = Field(..., description="Analysis configuration")
    
    # Results
    motifs: List[MotifRecord] = Field(default_factory=list, description="Detected motifs")
    
    # Summary statistics
    total_motifs: int = Field(0, ge=0, description="Total motifs detected")
    class_distribution: Dict[str, int] = Field(default_factory=dict, description="Motifs per class")
    coverage_statistics: Dict[str, float] = Field(default_factory=dict, description="Coverage statistics")
    
    # Performance metrics
    processing_time_seconds: Optional[float] = Field(None, ge=0.0, description="Processing time")
    memory_usage_mb: Optional[float] = Field(None, ge=0.0, description="Peak memory usage")
    
    @validator('total_motifs', always=True)
    def validate_motif_count(cls, v, values):
        """Ensure motif count matches actual motifs."""
        if 'motifs' in values:
            actual_count = len(values['motifs'])
            if v != actual_count:
                return actual_count
        return v

class ExportConfig(BaseModel):
    """Configuration for exporting results."""
    model_config = ConfigDict(str_strip_whitespace=True)
    
    format: Literal["bed", "csv", "parquet", "gff", "json"] = Field("bed", description="Export format")
    
    # BED format options
    include_score: bool = Field(True, description="Include score in BED output")
    score_type: Literal["raw", "normalized"] = Field("normalized", description="Score type for BED")
    color_by_class: bool = Field(True, description="Color motifs by class")
    
    # CSV/Parquet options
    include_sequence: bool = Field(True, description="Include sequence in tabular output")
    precision: int = Field(3, ge=1, le=10, description="Decimal precision for scores")
    
    # File options
    compress: bool = Field(False, description="Compress output files")
    output_directory: Optional[str] = Field(None, description="Output directory path")

class BatchAnalysisRequest(BaseModel):
    """Request for batch analysis of multiple sequences."""
    model_config = ConfigDict(str_strip_whitespace=True)
    
    sequences: Dict[str, str] = Field(..., description="Sequences to analyze")
    config: AnalysisConfig = Field(default_factory=AnalysisConfig, description="Analysis configuration")
    export_config: Optional[ExportConfig] = Field(None, description="Export configuration")
    
    # Batch options
    parallel_sequences: bool = Field(True, description="Process sequences in parallel")
    output_prefix: str = Field("nbdfinder_batch", description="Output file prefix")

class PerformanceMetrics(BaseModel):
    """Performance metrics for analysis."""
    model_config = ConfigDict(str_strip_whitespace=True)
    
    sequence_length: int = Field(..., ge=1, description="Sequence length analyzed")
    processing_time_seconds: float = Field(..., ge=0.0, description="Total processing time")
    memory_peak_mb: float = Field(..., ge=0.0, description="Peak memory usage")
    
    # Throughput metrics
    bases_per_second: float = Field(..., ge=0.0, description="Processing speed")
    motifs_per_second: float = Field(..., ge=0.0, description="Motif detection rate")
    
    # Algorithm breakdown
    pattern_matching_time: Optional[float] = Field(None, ge=0.0, description="Pattern matching time")
    scoring_time: Optional[float] = Field(None, ge=0.0, description="Scoring calculation time")
    postprocessing_time: Optional[float] = Field(None, ge=0.0, description="Post-processing time")

# Factory functions for creating common model instances

def create_motif_record(motif_class: str, subclass: str, start: int, end: int, 
                       sequence: str, score: float = 0.0, **kwargs) -> MotifRecord:
    """Factory function for creating motif records."""
    return MotifRecord(
        motif_class=MotifClass(motif_class),
        subclass=subclass,
        start=start,
        end=end,
        sequence=sequence,
        score=score,
        **kwargs
    )

def create_sequence_info(name: str, sequence: str, **kwargs) -> SequenceInfo:
    """Factory function for creating sequence info."""
    length = len(sequence)
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    gc_content = 100 * gc_count / length if length > 0 else 0.0
    
    return SequenceInfo(
        name=name,
        length=length,
        gc_content=gc_content,
        **kwargs
    )

__all__ = [
    # Enums
    'MotifClass',
    'Strand',
    
    # Models
    'MotifRecord',
    'SequenceInfo',
    'AnalysisConfig',
    'AnalysisResults',
    'ExportConfig',
    'BatchAnalysisRequest',
    'PerformanceMetrics',
    
    # Factories
    'create_motif_record',
    'create_sequence_info'
]