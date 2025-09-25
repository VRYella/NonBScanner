#!/usr/bin/env python3
"""
NBDScanner REST API
==================

Consolidated REST API for Non-B DNA Motif Detection
Supports all 11 major classes and 22+ subclasses using the consolidated NBDScanner engine

Endpoints:
- GET  /api/v1/health - Health check
- GET  /api/v1/classes - List all supported motif classes  
- GET  /api/v1/patterns - Get pattern statistics
- POST /api/v1/analyze - Analyze DNA sequence for all motifs
- GET  /api/v1/motif-info - Get comprehensive motif information
- GET  /api/v1/stats - Get API usage statistics
"""

import os
import sys
import time
import logging
from typing import List, Dict, Any, Optional
from datetime import datetime

import uvicorn
from fastapi import FastAPI, HTTPException, Request, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from contextlib import asynccontextmanager

# Add current directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import consolidated NBDScanner modules
from nbdscanner import (
    analyze_sequence, analyze_multiple_sequences, 
    get_motif_classification_info
)
from utils_consolidated import (
    parse_fasta, validate_sequence, get_basic_stats,
    export_to_json, export_to_csv, export_to_bed
)
from motif_patterns import (
    PatternRegistry, get_pattern_statistics
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# API Statistics storage
api_stats = {
    "total_requests": 0,
    "analyze_requests": 0,
    "class_requests": 0,
    "start_time": datetime.utcnow(),
    "motifs_detected": 0,
    "sequences_analyzed": 0
}

@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan management"""
    logger.info("NBDScanner API starting up...")
    yield
    logger.info("NBDScanner API shutting down...")

# Initialize FastAPI app
app = FastAPI(
    title="NBDScanner API", 
    description="Consolidated REST API for Non-B DNA Motif Detection with 11 classes and 22+ subclasses",
    version="2024.1",
    docs_url="/docs",
    redoc_url="/redoc", 
    lifespan=lifespan
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Request/Response Models
class SequenceRequest(BaseModel):
    sequence: str = Field(..., description="DNA sequence to analyze", min_length=10)
    sequence_name: str = Field(default="sequence", description="Name identifier for the sequence")
    detailed: bool = Field(default=True, description="Include detailed analysis metadata")

class MotifResult(BaseModel):
    Class: str
    Subclass: str
    Start: int
    End: int
    Length: int
    Sequence: str
    Score: float
    Strand: str
    Method: str

class AnalysisResponse(BaseModel):
    sequence_name: str
    sequence_length: int
    total_motifs: int
    classes_detected: int
    subclasses_detected: int
    analysis_time: float
    motifs: List[MotifResult]
    length_limits: Dict[str, int]

class MotifResult(BaseModel):
    sequence_name: str
    total_motifs: int
    classes_detected: List[str]
    subclasses_detected: List[str]
    analysis_time: float
    motifs: List[Dict[str, Any]]

class HealthResponse(BaseModel):
    status: str
    version: str
    timestamp: str
    available_detectors: List[str]

class StatsResponse(BaseModel):
    total_requests: int
    analyze_requests: int
    class_requests: int
    motifs_detected: int
    sequences_analyzed: int
    uptime_seconds: float
    start_time: str

# Middleware for request tracking
@app.middleware("http")
async def track_requests(request: Request, call_next):
    global api_stats
    api_stats["total_requests"] += 1
    
    start_time = time.time()
    response = await call_next(request)
    process_time = time.time() - start_time
    
    response.headers["X-Process-Time"] = str(process_time)
    return response

# API Endpoints

@app.get("/api/v1/health", response_model=HealthResponse)
async def health_check():
    """Health check endpoint"""
    return HealthResponse(
        status="healthy",
        version="1.0.0",
        timestamp=datetime.utcnow().isoformat(),
        available_detectors=get_available_detectors()
    )

@app.get("/api/v1/classes", response_model=List[MotifClassInfo])
async def get_motif_classes():
    """Get information about all supported motif classes"""
    global api_stats
    api_stats["class_requests"] += 1
    
    classes_info = []
    
    # Define comprehensive class information including subclasses
    class_definitions = {
        "Curved_DNA": {
            "name": "Curved DNA",
            "description": "A-tract mediated DNA curvature",
            "subclasses": ["Global Curvature", "Local Curvature", "Phased A-tracts"]
        },
        "Slipped_DNA": {
            "name": "Slipped DNA", 
            "description": "Direct/tandem repeats forming slipped structures",
            "subclasses": ["Direct Repeat", "STR (Short Tandem Repeat)", "Slipped DNA [STR]"]
        },
        "Cruciform": {
            "name": "Cruciform DNA",
            "description": "Inverted repeats forming four-way junctions", 
            "subclasses": ["Palindromic Inverted Repeat", "Cruciform Hairpin"]
        },
        "R-Loop": {
            "name": "R-Loop",
            "description": "RNA-DNA hybrids with displaced ssDNA",
            "subclasses": ["R-Loop Formation Site", "Stable R-Loop"]
        },
        "Triplex": {
            "name": "Triplex DNA",
            "description": "Three-stranded DNA structures",
            "subclasses": ["Mirror_Repeat_Triplex", "G-Triplex intermediate", "Sticky DNA"]
        },
        "G-Quadruplex": {
            "name": "G-Quadruplex Family",
            "description": "G-quadruplex structures and variants",
            "subclasses": [
                "Canonical G4", "Bulged G4", "Relaxed G4", "Bipartite G4", 
                "Multimeric G4", "Imperfect G4", "G-Triplex"
            ]
        },
        "i-Motif": {
            "name": "i-Motif",
            "description": "C-rich structures complementary to G4",
            "subclasses": ["Canonical i-Motif", "Extended i-Motif", "AC-Motif"]
        },
        "Z-DNA": {
            "name": "Z-DNA",
            "description": "Left-handed double helix",
            "subclasses": ["Z-DNA", "eGZ (Extruded-G Z-DNA)", "GC-rich Z-DNA"]
        },
        "Hybrid": {
            "name": "Hybrid Motifs",
            "description": "Overlapping/composite motifs",
            "subclasses": [
                "G-Quadruplex_Triplex_DNA_Overlap", "Slipped_DNA_Z-DNA_Overlap",
                "Multi-class Overlap"
            ]
        },
        "Cluster": {
            "name": "Non-B DNA Clusters", 
            "description": "High-density motif regions",
            "subclasses": ["Motif Hotspot", "Dense Cluster", "Mixed Cluster"]
        }
    }
    
    for class_id, class_info in class_definitions.items():
        # Get scoring method info
        scoring_info = SCORING_METHODS.get(class_id, SCORING_METHODS.get("Curved_DNA", {}))
        scoring_method = scoring_info.get("method", "default")
        
        # Get length limits
        length_limits = MOTIF_LENGTH_LIMITS.get(class_id, {"S_min": 10, "S_max": 200})
        
        classes_info.append(MotifClassInfo(
            class_id=class_id,
            class_name=class_info["name"],
            description=class_info["description"],
            subclasses=class_info["subclasses"],
            scoring_method=scoring_method,
            length_limits=length_limits
        ))
    
    return classes_info

@app.get("/api/v1/classes/{class_id}", response_model=MotifClassInfo)
async def get_motif_class_details(class_id: str):
    """Get detailed information about a specific motif class"""
    global api_stats
    api_stats["class_requests"] += 1
    
    # Get all classes and find the requested one
    all_classes = await get_motif_classes()
    
    for class_info in all_classes:
        if class_info.class_id == class_id:
            return class_info
    
    raise HTTPException(
        status_code=404,
        detail=f"Motif class '{class_id}' not found"
    )

@app.post("/api/v1/analyze", response_model=AnalysisResponse)
async def analyze_dna_sequence(request: SequenceRequest, background_tasks: BackgroundTasks):
    """Analyze DNA sequence using consolidated NBDScanner"""
    global api_stats
    api_stats["analyze_requests"] += 1
    api_stats["sequences_analyzed"] += 1
    
    start_time = time.time()
    
    try:
        # Validate sequence using consolidated utilities
        sequence = request.sequence.upper().strip()
        is_valid, error_msg = validate_sequence(sequence)
        if not is_valid:
            raise HTTPException(status_code=400, detail=error_msg)
        
        # Run consolidated NBDScanner analysis
        motifs = analyze_sequence(
            sequence=sequence,
            sequence_name=request.sequence_name,
            detailed=request.detailed
        )
        
        # Calculate analysis statistics
        total_motifs = len(motifs)
        classes_detected = len(set(m.get('Class', 'Unknown') for m in motifs))
        subclasses_detected = len(set(m.get('Subclass', 'Unknown') for m in motifs))
        analysis_time = time.time() - start_time
        
        # Update API statistics
        api_stats["motifs_detected"] += total_motifs
        
        # Convert to response model
        motif_results = [MotifResult(**motif) for motif in motifs]
        
        return AnalysisResponse(
            sequence_name=request.sequence_name,
            sequence_length=len(sequence),
            total_motifs=total_motifs,
            classes_detected=classes_detected,
            subclasses_detected=subclasses_detected,
            analysis_time=round(analysis_time, 3),
            motifs=motif_results
        )
        
    except Exception as e:
        logger.error(f"Analysis error: {str(e)}")
        raise HTTPException(
            status_code=500,
            detail=f"Analysis failed: {str(e)}"
        )

@app.post("/api/v1/analyze/{class_id}")
async def analyze_sequence_class(class_id: str, request: SequenceRequest):
    """Analyze DNA sequence for specific motif class"""
    global api_stats
    api_stats["analyze_requests"] += 1
    
    start_time = time.time()
    
    try:
        # Validate sequence
        sequence = request.sequence.upper().strip()
        if not all(c in 'ATGCN' for c in sequence):
            raise HTTPException(
                status_code=400,
                detail="Invalid DNA sequence. Only ATGCN characters allowed."
            )
        
        # Check if class exists
        if class_id not in DETECTOR_REGISTRY:
            available_classes = list(DETECTOR_REGISTRY.keys())
            raise HTTPException(
                status_code=404,
                detail=f"Class '{class_id}' not available. Available: {available_classes}"
            )
        
        # Run class-specific analysis
        motifs = detect_motif_class(
            sequence=sequence,
            motif_class=class_id,
            sequence_name=request.sequence_name
        )
        
        analysis_time = time.time() - start_time
        
        return {
            "sequence_name": request.sequence_name,
            "class_analyzed": class_id,
            "total_motifs": len(motifs),
            "analysis_time": round(analysis_time, 3),
            "motifs": motifs
        }
        
    except Exception as e:
        logger.error(f"Class analysis error: {str(e)}")
        raise HTTPException(
            status_code=500,
            detail=f"Class analysis failed: {str(e)}"
        )

@app.get("/api/v1/stats", response_model=StatsResponse)
async def get_api_statistics():
    """Get API usage statistics"""
    uptime = (datetime.utcnow() - api_stats["start_time"]).total_seconds()
    
    return StatsResponse(
        total_requests=api_stats["total_requests"],
        analyze_requests=api_stats["analyze_requests"],
        class_requests=api_stats["class_requests"],
        motifs_detected=api_stats["motifs_detected"],
        sequences_analyzed=api_stats["sequences_analyzed"],
        uptime_seconds=round(uptime, 2),
        start_time=api_stats["start_time"].isoformat()
    )

@app.get("/api/v1/motif-info")
async def get_comprehensive_motif_info():
    """Get comprehensive information about all 11 motif classes and 22+ subclasses"""
    
    # Use consolidated classification info
    classification_info = get_motif_classification_info()
    pattern_stats = get_pattern_statistics()
    
    comprehensive_info = {
        "overview": {
            "version": classification_info["version"],
            "total_classes": classification_info["total_classes"],
            "total_subclasses": classification_info["total_subclasses"],
            "total_patterns": pattern_stats["total_patterns"],
            "description": "NBDScanner detects Non-B DNA structures across 11 major classes with 22+ specialized subclasses using 52+ regex patterns"
        },
        "classification": classification_info["classification"],
        "pattern_statistics": pattern_stats,
        "technology": {
            "engine": "NBDScanner Consolidated Engine",
            "pattern_matching": "Hyperscan accelerated (when available)",
            "scoring_algorithms": pattern_stats["scoring_methods"],
            "performance": "40x+ speed improvement over traditional methods",
            "output_format": "Standardized genomic coordinates with comprehensive metadata"
        },
        "references": [
            "Bedrat et al. (2016) Nucleic Acids Research - G4Hunter algorithm", 
            "Ho et al. (1986) Nucleic Acids Research - Z-DNA scoring",
            "Olson et al. (1998) PNAS - DNA curvature models",
            "Frank-Kamenetskii & Mirkin (1995) - Triplex DNA",
            "Zeraati et al. (2018) Nature Chemistry - i-motif structures"
        ]
    }
    
    return comprehensive_info
                "biological_relevance": "Recombination hotspots, gene regulation, genomic instability"
            },
            "Class_04_R_Loop": {
                "name": "R-Loop",
                "description": "RNA-DNA hybrids with displaced single-stranded DNA",
                "subclasses": [
                    {"id": "4.1", "name": "R-Loop Formation Site", "description": "Regions prone to R-loop formation"}
                ],
                "scoring_method": "RLFS model considering G-richness and transcriptional activity",
                "biological_relevance": "Transcription regulation, genome instability, DNA damage"
            },
            "Class_05_Triplex": {
                "name": "Triplex DNA",
                "description": "Three-stranded DNA structures",
                "subclasses": [
                    {"id": "5.1", "name": "Mirror Repeat Triplex", "description": "Homopurine/homopyrimidine mirror repeats"},
                    {"id": "5.2", "name": "Sticky DNA", "description": "Extended GAA/TTC repeats forming triplex structures"}
                ],
                "scoring_method": "Triplex potential based on homogeneity and pH stability",
                "biological_relevance": "Gene regulation, genetic diseases, therapeutic targets"
            },
            "Class_06_G4_Family": {
                "name": "G-Quadruplex Family",
                "description": "G-rich sequences forming four-stranded structures",
                "subclasses": [
                    {"id": "6.1", "name": "Canonical G4", "description": "Classic G3+N1-7G3+N1-7G3+N1-7G3+ pattern"},
                    {"id": "6.2", "name": "Bulged G4", "description": "G4 structures with bulge loops"},
                    {"id": "6.3", "name": "Relaxed G4", "description": "G2+ patterns with relaxed criteria"},
                    {"id": "6.4", "name": "Bipartite G4", "description": "Two G4-forming regions with long spacers"},
                    {"id": "6.5", "name": "Multimeric G4", "description": "Multiple G4 units in tandem"},
                    {"id": "6.6", "name": "Imperfect G4", "description": "G4-like structures with interruptions"},
                    {"id": "6.7", "name": "G-Triplex", "description": "Three G-tracts forming intermediate structures"}
                ],
                "scoring_method": "G4Hunter algorithm with class-specific thresholds",
                "biological_relevance": "Telomere biology, gene regulation, cancer therapy targets"
            },
            "Class_07_i_Motif": {
                "name": "i-Motif",
                "description": "C-rich structures complementary to G-quadruplexes",
                "subclasses": [
                    {"id": "7.1", "name": "Canonical i-Motif", "description": "C4+N1-7C4+N1-7C4+N1-7C4+ patterns"},
                    {"id": "7.2", "name": "Extended i-Motif", "description": "Longer C-rich regions"},
                    {"id": "7.3", "name": "AC-Motif", "description": "Alternating A-rich/C-rich consensus regions"}
                ],
                "scoring_method": "Adapted G4Hunter for C-richness with pH considerations",
                "biological_relevance": "Gene regulation, pH sensing, G4 regulation"
            },
            "Class_08_Z_DNA": {
                "name": "Z-DNA",
                "description": "Left-handed double helix structures",
                "subclasses": [
                    {"id": "8.1", "name": "Z-DNA", "description": "Classic alternating purine-pyrimidine Z-DNA"},
                    {"id": "8.2", "name": "eGZ (Extruded-G Z-DNA)", "description": "CGG repeats with Z-DNA potential"}
                ],
                "scoring_method": "Z-seeker algorithm with dinucleotide scoring",
                "biological_relevance": "Gene regulation, chromatin structure, genetic diseases"
            },
            "Class_09_Hybrid": {
                "name": "Hybrid Motifs",
                "description": "Overlapping or composite motif structures",
                "subclasses": [
                    {"id": "9.1", "name": "Multi-class Overlap", "description": "Regions where multiple motif classes overlap"},
                    {"id": "9.2", "name": "Composite Structures", "description": "Complex arrangements of multiple motif types"}
                ],
                "scoring_method": "Overlap diversity and structural complexity scoring",
                "biological_relevance": "Complex regulatory regions, structural hotspots"
            },
            "Class_10_Cluster": {
                "name": "Non-B DNA Clusters",
                "description": "High-density regions with multiple motifs",
                "subclasses": [
                    {"id": "10.1", "name": "Motif Hotspot", "description": "Regions with high motif density"},
                    {"id": "10.2", "name": "Mixed Cluster", "description": "Clusters with diverse motif types"}
                ],
                "scoring_method": "Density and diversity-based clustering analysis",
                "biological_relevance": "Regulatory hotspots, structural domains, genome organization"
            }
        },
        "technical_details": {
            "detection_method": "Hyperscan-accelerated pattern matching with scientific scoring",
            "scoring_normalization": "0-1 scale for cross-class comparison",
            "performance": "40x+ speed improvement over traditional methods",
            "output_format": "Standardized genomic coordinates with comprehensive metadata"
        },
        "references": [
            "Bedrat et al. (2016) Nucleic Acids Research - G4Hunter algorithm",
            "Ho et al. (1986) Nucleic Acids Research - Z-DNA scoring",
            "Olson et al. (1998) PNAS - DNA curvature models",
            "Crothers et al. (1992) Methods Enzymol - Phasing assays"
        ]
    }
    
    return comprehensive_info

# Development server
if __name__ == "__main__":
    uvicorn.run(
        "api:app",
        host="0.0.0.0",
        port=8000,
        reload=True,
        log_level="info"
    )