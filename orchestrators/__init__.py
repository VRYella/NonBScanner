"""
NBDFinder Orchestrators
=====================

High-level orchestration modules for coordinating motif detection across
multiple detectors and managing large-scale analysis workflows.
"""

# Import main orchestration functions
try:
    from .all_motifs import all_motifs as detect_all_motifs
except ImportError:
    detect_all_motifs = None

try:
    from .all_motifs import analyze_sequence_comprehensive
except ImportError:
    analyze_sequence_comprehensive = None

from .stream_orchestrator import StreamingOrchestrator, BatchStreamingOrchestrator

__all__ = [
    'detect_all_motifs',
    'analyze_sequence_comprehensive',
    'StreamingOrchestrator',
    'BatchStreamingOrchestrator'
]