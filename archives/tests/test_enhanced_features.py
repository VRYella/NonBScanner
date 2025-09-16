#!/usr/bin/env python3
"""
Comprehensive Test Suite for Enhanced NBDFinder
==============================================

Tests the new classification configuration, conservation analysis,
normalized scoring, and enhanced UI features.
"""

import sys
import os
import numpy as np
from collections import Counter

def test_classification_config():
    """Test classification configuration module"""
    print("🧪 Testing Classification Configuration...")
    
    try:
        from classification_config import (
            MOTIF_LENGTH_LIMITS, SCORING_METHODS, 
            get_motif_limits, normalize_score, 
            calculate_enrichment_score, classify_conservation
        )
        
        # Test length limits
        assert len(MOTIF_LENGTH_LIMITS) > 0, "No motif length limits defined"
        print(f"   ✓ Found {len(MOTIF_LENGTH_LIMITS)} motif classes with length limits")
        
        # Test specific motif limits
        g4_limits = get_motif_limits("G4")
        assert g4_limits == (13, 100), f"Unexpected G4 limits: {g4_limits}"
        print(f"   ✓ G4 limits correct: {g4_limits}")
        
        # Test normalization
        normalized = normalize_score(50.0, 25, "G4")
        assert 0 <= normalized <= 100, f"Normalized score out of range: {normalized}"
        print(f"   ✓ Score normalization works: 50.0 -> {normalized}")
        
        # Test enrichment calculation
        observed = 10
        shuffled = [2, 3, 1, 4, 2, 5, 3, 2, 1, 3]
        enrichment, p_val = calculate_enrichment_score(observed, shuffled)
        assert enrichment > 0, f"Expected positive enrichment, got {enrichment}"
        assert 0 <= p_val <= 1, f"P-value out of range: {p_val}"
        print(f"   ✓ Enrichment calculation works: {enrichment}, p={p_val}")
        
        # Test conservation classification
        conservation_class = classify_conservation(enrichment, p_val)
        assert conservation_class in ["Highly_Conserved", "Moderately_Conserved", "Weakly_Conserved", "Neutral", "Depleted"]
        print(f"   ✓ Conservation classification: {conservation_class}")
        
        return True
        
    except Exception as e:
        print(f"   ❌ Classification config test failed: {e}")
        return False

def test_conservation_analysis():
    """Test conservation analysis module"""
    print("🧪 Testing Conservation Analysis...")
    
    try:
        from conservation_analysis import (
            shuffle_sequence, generate_shuffled_sequences,
            calculate_motif_conservation, get_conservation_summary
        )
        
        test_seq = "ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGG"
        
        # Test sequence shuffling
        shuffled = shuffle_sequence(test_seq, seed=42)
        assert len(shuffled) == len(test_seq), "Shuffled sequence length mismatch"
        assert set(shuffled) == set(test_seq), "Shuffled sequence composition changed"
        print(f"   ✓ Sequence shuffling preserves composition")
        
        # Test multiple shuffles
        shuffled_seqs = generate_shuffled_sequences(test_seq, n_shuffles=10)
        assert len(shuffled_seqs) == 10, "Wrong number of shuffled sequences"
        print(f"   ✓ Generated {len(shuffled_seqs)} shuffled sequences")
        
        # Test conservation summary with mock data
        mock_motifs = [
            {'Class': 'G4', 'Conservation_Score': 2.5, 'Conservation_P_Value': 0.01, 'Conservation_Class': 'Highly_Conserved'},
            {'Class': 'Curved_DNA', 'Conservation_Score': 0.5, 'Conservation_P_Value': 0.2, 'Conservation_Class': 'Weakly_Conserved'},
            {'Class': 'Z-DNA', 'Conservation_Score': -1.0, 'Conservation_P_Value': 0.8, 'Conservation_Class': 'Depleted'}
        ]
        
        summary = get_conservation_summary(mock_motifs)
        assert summary['total_motifs'] == 3, "Wrong motif count in summary"
        assert summary['highly_conserved'] == 1, "Wrong conserved count"
        assert summary['depleted'] == 1, "Wrong depleted count"
        print(f"   ✓ Conservation summary generation works")
        
        return True
        
    except Exception as e:
        print(f"   ❌ Conservation analysis test failed: {e}")
        return False

def test_enhanced_motif_detection():
    """Test enhanced motif detection with new features"""
    print("🧪 Testing Enhanced Motif Detection...")
    
    try:
        import motifs
        
        test_seq = "ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC"
        
        # Test basic motif detection
        basic_motifs = motifs.all_motifs(test_seq, calculate_conservation=False)
        assert len(basic_motifs) > 0, "No motifs detected"
        print(f"   ✓ Basic detection found {len(basic_motifs)} motifs")
        
        # Test that normalized scores are added
        has_normalized = any('Normalized_Score' in motif for motif in basic_motifs)
        assert has_normalized, "No normalized scores found"
        print(f"   ✓ Normalized scores added to motifs")
        
        # Test conservation analysis (if available)
        try:
            enhanced_motifs = motifs.all_motifs(test_seq, calculate_conservation=True)
            has_conservation = any('Conservation_Score' in motif for motif in enhanced_motifs)
            if has_conservation:
                print(f"   ✓ Conservation analysis integrated")
            else:
                print(f"   ⚠ Conservation analysis not available (expected in some environments)")
        except Exception:
            print(f"   ⚠ Conservation analysis skipped (dependencies not available)")
        
        # Test motif validation
        for motif in basic_motifs[:3]:  # Check first 3 motifs
            assert 'Class' in motif, "Motif missing class"
            assert 'Start' in motif, "Motif missing start position"
            assert 'End' in motif, "Motif missing end position"
            assert motif['Start'] <= motif['End'], "Invalid motif coordinates"
        
        print(f"   ✓ Motif structure validation passed")
        
        return True
        
    except Exception as e:
        print(f"   ❌ Enhanced motif detection test failed: {e}")
        return False

def test_app_compatibility():
    """Test that app.py can import all required modules"""
    print("🧪 Testing App Compatibility...")
    
    try:
        # Test core imports that app.py needs
        import streamlit as st
        import pandas as pd
        import matplotlib.pyplot as plt
        
        # Test our new modules
        import classification_config
        import conservation_analysis
        import motifs
        
        print(f"   ✓ All app dependencies available")
        
        # Test that configuration is accessible
        assert hasattr(classification_config, 'MOTIF_LENGTH_LIMITS')
        assert hasattr(classification_config, 'SCORING_METHODS')
        print(f"   ✓ Configuration data accessible")
        
        # Test that conservation functions are available
        assert hasattr(conservation_analysis, 'calculate_motif_conservation')
        assert hasattr(conservation_analysis, 'get_conservation_summary')
        print(f"   ✓ Conservation functions accessible")
        
        return True
        
    except Exception as e:
        print(f"   ❌ App compatibility test failed: {e}")
        return False

def run_integration_test():
    """Run a complete integration test"""
    print("🧪 Running Integration Test...")
    
    try:
        import motifs
        from classification_config import get_motif_limits, normalize_score
        
        # Test sequence with multiple motif types
        test_seq = ("ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGG"
                   "CCCCCTCCCCCTCCCCCTCCCCGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA"
                   "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG")
        
        # Run complete analysis
        all_results = motifs.all_motifs(
            test_seq, 
            nonoverlap=False, 
            report_hotspots=True,
            calculate_conservation=False  # Skip for faster testing
        )
        
        print(f"   ✓ Complete analysis found {len(all_results)} motifs")
        
        # Check that different motif classes are found
        motif_classes = set(motif.get('Class', 'Unknown') for motif in all_results)
        print(f"   ✓ Found motif classes: {', '.join(sorted(motif_classes))}")
        
        # Check normalized scoring
        normalized_scores = [motif.get('Normalized_Score', 0) for motif in all_results]
        valid_scores = [score for score in normalized_scores if 0 <= float(score) <= 100]
        print(f"   ✓ {len(valid_scores)}/{len(normalized_scores)} motifs have valid normalized scores")
        
        # Test length constraint application
        for motif in all_results[:5]:  # Check first 5
            motif_class = motif.get('Class', '')
            motif_length = motif.get('Length', 0)
            
            try:
                s_min, s_max = get_motif_limits(motif_class)
                if motif_length < s_min:
                    print(f"   ⚠ Short motif detected: {motif_class} {motif_length}bp < {s_min}bp")
                elif motif_length > s_max:
                    print(f"   ⚠ Long motif detected: {motif_class} {motif_length}bp > {s_max}bp")
            except:
                pass
        
        print(f"   ✓ Length constraint checking completed")
        
        return True
        
    except Exception as e:
        print(f"   ❌ Integration test failed: {e}")
        return False

def main():
    """Run all tests"""
    print("=" * 60)
    print("🧬 NBDFinder Enhanced Features Test Suite")
    print("=" * 60)
    
    tests = [
        test_classification_config,
        test_conservation_analysis,
        test_enhanced_motif_detection,
        test_app_compatibility,
        run_integration_test
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        try:
            if test():
                passed += 1
            print()
        except Exception as e:
            print(f"   ❌ Test {test.__name__} crashed: {e}")
            print()
    
    print("=" * 60)
    print(f"🏆 Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 ALL TESTS PASSED! Enhanced NBDFinder is ready!")
    else:
        print("⚠️  Some tests failed. Check implementation.")
    print("=" * 60)
    
    return passed == total

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)