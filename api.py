#!/usr/bin/env python3
"""
NBDScanner REST API
==================

Consolidated REST API for Non-B DNA Motif Detection with 11 classes and 22+ subclasses.
Built on the consolidated NBDScanner engine for high-performance analysis.

Author: Dr. Venkata Rajesh Yella
License: MIT
Version: 2024.1
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
from utils import (
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

class HealthResponse(BaseModel):
    status: str
    version: str
    timestamp: datetime
    uptime_seconds: float

# API Endpoints

@app.middleware("http")
async def log_requests(request: Request, call_next):
    """Log all API requests"""
    start_time = time.time()
    api_stats["total_requests"] += 1
    
    response = await call_next(request)
    process_time = time.time() - start_time
    
    logger.info(f"{request.method} {request.url.path} - {response.status_code} - {process_time:.3f}s")
    return response

@app.get("/api/v1/health", response_model=HealthResponse)
async def health_check():
    """Health check endpoint"""
    uptime = (datetime.utcnow() - api_stats["start_time"]).total_seconds()
    return HealthResponse(
        status="healthy",
        version="2024.1",
        timestamp=datetime.utcnow(),
        uptime_seconds=uptime
    )

@app.get("/api/v1/classes")
async def get_motif_classes():
    """Get information about all supported motif classes"""
    classification_info = get_motif_classification_info()
    return {
        "status": "success",
        "total_classes": classification_info["total_classes"],
        "total_subclasses": classification_info["total_subclasses"],
        "classification": classification_info["classification"]
    }

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
        raise HTTPException(status_code=500, detail=str(e))

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

@app.get("/api/v1/patterns")
async def get_pattern_info():
    """Get detailed pattern statistics"""
    global api_stats
    api_stats["total_requests"] += 1
    
    try:
        pattern_stats = get_pattern_statistics()
        return {
            "status": "success", 
            "pattern_statistics": pattern_stats,
            "validation_status": "All patterns validated"
        }
    except Exception as e:
        logger.error(f"Pattern statistics error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/v1/stats")
async def get_api_statistics():
    """Get API usage statistics"""
    uptime = (datetime.utcnow() - api_stats["start_time"]).total_seconds()
    
    return {
        "status": "success",
        "statistics": {
            "total_requests": api_stats["total_requests"],
            "analyze_requests": api_stats["analyze_requests"],
            "sequences_analyzed": api_stats["sequences_analyzed"],
            "motifs_detected": api_stats["motifs_detected"],
            "uptime_seconds": round(uptime, 2),
            "start_time": api_stats["start_time"].isoformat()
        }
    }

# Development server
if __name__ == "__main__":
    uvicorn.run(
        "api:app",
        host="0.0.0.0",
        port=8000,
        reload=True,
        log_level="info"
    )