"""
NBDFinder Command Line Interface
===============================

Command-line tools for NBDFinder including:
- Main analysis CLI
- FASTA slicing utilities
"""

from .main import main as cli_main
from .slice_fasta import main as slice_main

__all__ = [
    'cli_main',
    'slice_main'
]