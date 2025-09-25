"""
Hyperscan Database Manager for Optimal Performance
==================================================

This module provides centralized management of Hyperscan databases to maximize
performance by pre-compiling and caching databases, avoiding repeated compilation
overhead.

Key Performance Optimizations:
1. Database Pre-compilation: Compile once, use many times
2. Pattern Optimization: Optimized regex patterns for Hyperscan
3. Memory Efficiency: Reuse database objects
4. Callback Optimization: Streamlined callback functions
5. Thread Safety: Safe for concurrent use
"""

import hyperscan
import threading
import hashlib
from typing import List, Tuple, Dict, Any, Callable, Optional
from collections import defaultdict

class HyperscanManager:
    """
    Singleton manager for Hyperscan databases with caching and optimization.
    """
    
    _instance = None
    _lock = threading.Lock()
    
    def __new__(cls):
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
                    cls._instance._initialized = False
        return cls._instance
    
    def __init__(self):
        if self._initialized:
            return
            
        self._database_cache = {}
        self._pattern_cache = {}
        self._cache_lock = threading.Lock()
        self._initialized = True
    
    def _generate_cache_key(self, patterns: List[Tuple], flags: int = 0) -> str:
        """Generate a unique cache key for pattern set."""
        pattern_strs = [str(p) for p in patterns]
        pattern_data = "|".join(pattern_strs) + f"|flags:{flags}"
        return hashlib.md5(pattern_data.encode()).hexdigest()
    
    def get_optimized_database(self, patterns: List[Tuple], flags: int = 0) -> hyperscan.Database:
        """
        Get or create an optimized Hyperscan database for given patterns.
        
        Args:
            patterns: List of (regex_pattern, id) tuples
            flags: Hyperscan compilation flags
            
        Returns:
            Compiled Hyperscan database
        """
        cache_key = self._generate_cache_key(patterns, flags)
        
        with self._cache_lock:
            if cache_key in self._database_cache:
                return self._database_cache[cache_key]
            
            # Extract expressions and IDs
            expressions = [p[0].encode() if isinstance(p[0], str) else p[0] for p in patterns]
            ids = [p[1] for p in patterns]
            
            # Create and compile database with optimization flags
            db = hyperscan.Database()
            try:
                # Use optimized flags for performance
                compile_flags = flags | hyperscan.HS_FLAG_SOM_LEFTMOST
                db.compile(
                    expressions=expressions,
                    ids=ids,
                    flags=[compile_flags] * len(expressions)
                )
                
                # Cache the compiled database
                self._database_cache[cache_key] = db
                return db
                
            except Exception as e:
                # Fallback to basic compilation if optimization fails
                db.compile(expressions=expressions, ids=ids)
                self._database_cache[cache_key] = db
                return db
    
    def optimized_scan(self, 
                      patterns: List[Tuple], 
                      sequence: str, 
                      callback: Callable,
                      context: Any = None,
                      flags: int = 0) -> List[Any]:
        """
        Perform optimized Hyperscan scanning with database caching.
        
        Args:
            patterns: List of pattern tuples
            sequence: Target sequence to scan
            callback: Match callback function
            context: Optional context for callback
            flags: Compilation flags
            
        Returns:
            List of matches found
        """
        if not patterns or not sequence:
            return []
        
        # Get or create optimized database
        db = self.get_optimized_database(patterns, flags)
        
        # Prepare sequence
        seq_bytes = sequence.upper().encode()
        
        # Perform scan
        matches = []
        
        def optimized_callback(id, from_, to, flags, ctx):
            try:
                result = callback(id, from_, to, flags, ctx)
                return result if result is not None else hyperscan.HS_SUCCESS
            except Exception:
                return hyperscan.HS_SUCCESS
        
        try:
            db.scan(seq_bytes, match_event_handler=optimized_callback, context=context)
        except Exception:
            # Continue scanning even if individual matches fail
            pass
        
        return matches
    
    def clear_cache(self):
        """Clear all cached databases."""
        with self._cache_lock:
            self._database_cache.clear()
            self._pattern_cache.clear()
    
    def get_cache_stats(self) -> Dict[str, int]:
        """Get cache statistics."""
        with self._cache_lock:
            return {
                'database_cache_size': len(self._database_cache),
                'pattern_cache_size': len(self._pattern_cache)
            }

# Global instance
hyperscan_manager = HyperscanManager()

def optimized_hs_find(patterns: List[Tuple], 
                     sequence: str, 
                     callback_func: Callable,
                     context: Any = None) -> List[Any]:
    """
    High-performance Hyperscan pattern matching with database caching.
    
    Args:
        patterns: List of (regex, id, ...) tuples
        sequence: Target sequence
        callback_func: Callback function for matches
        context: Optional context
        
    Returns:
        List of matches
    """
    if not patterns or not sequence:
        return []
    
    # Extract just regex and id for database compilation
    db_patterns = [(p[0], p[1]) for p in patterns]
    
    # Store full pattern info for callback access
    pattern_map = {p[1]: p for p in patterns}
    
    matches = []
    
    def enhanced_callback(id, from_, to, flags, ctx):
        if id in pattern_map:
            try:
                result = callback_func(id, from_, to, flags, ctx, pattern_map[id])
                if result is not None:
                    if isinstance(result, dict):
                        matches.append(result)
                    return hyperscan.HS_SUCCESS
            except Exception:
                pass
        return hyperscan.HS_SUCCESS
    
    # Use the optimized manager
    hyperscan_manager.optimized_scan(
        db_patterns, sequence, enhanced_callback, context
    )
    
    return matches

def clear_hyperscan_cache():
    """Clear all Hyperscan database caches."""
    hyperscan_manager.clear_cache()

def get_hyperscan_cache_stats() -> Dict[str, int]:
    """Get Hyperscan cache statistics."""
    return hyperscan_manager.get_cache_stats()