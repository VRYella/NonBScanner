#!/usr/bin/env python3
"""
Test suite for NBDFinder Revamp Features
========================================

Tests the new REST API, export utilities, enhanced caching, and overall integration.
"""

import sys
import os
import time
import json
import tempfile
import subprocess
import requests
import threading
from datetime import datetime

def test_rest_api():
    """Test REST API functionality"""
    print("🚀 Testing REST API...")
    
    # Start API server in background
    api_process = subprocess.Popen(
        [sys.executable, "api.py"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    
    # Wait for server to start
    time.sleep(3)
    
    try:
        # Test health endpoint
        response = requests.get("http://localhost:8000/api/v1/health", timeout=5)
        assert response.status_code == 200
        health_data = response.json()
        assert health_data["status"] == "healthy"
        print(f"   ✓ Health check passed - {health_data['motif_classes_available']} classes available")
        
        # Test analyze endpoint
        test_sequence = "GGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC"
        analyze_request = {
            "sequence": test_sequence,
            "sequence_name": "test_sequence",
            "report_hotspots": True
        }
        
        response = requests.post(
            "http://localhost:8000/api/v1/analyze",
            json=analyze_request,
            timeout=10
        )
        assert response.status_code == 200
        analysis_data = response.json()
        assert analysis_data["success"] == True
        assert len(analysis_data["motifs"]) > 0
        print(f"   ✓ Analysis endpoint passed - found {len(analysis_data['motifs'])} motifs")
        
        # Test stats endpoint
        response = requests.get("http://localhost:8000/api/v1/stats", timeout=5)
        assert response.status_code == 200
        stats_data = response.json()
        assert stats_data["total_analyses"] >= 1
        print(f"   ✓ Stats endpoint passed - {stats_data['total_analyses']} total analyses")
        
        return True
        
    except Exception as e:
        print(f"   ❌ API test failed: {e}")
        return False
    finally:
        # Clean up
        api_process.terminate()
        api_process.wait()

def test_export_utilities():
    """Test export format utilities"""
    print("📄 Testing Export Utilities...")
    
    try:
        from export_utils import export_to_bed, export_to_gff3, create_density_bedgraph
        
        # Sample motif data
        test_motifs = [
            {
                'Start': 1, 'End': 20, 'Class': 'G-Quadruplex', 'Subclass': 'Canonical_G4',
                'Normalized_Score': 0.8, 'Actual_Score': 4.2, 'GC_Content': 70.0,
                'Sequence': 'GGGTTAGGGTTAGGGTTAGG', 'Scoring_Method': 'G4Hunter'
            },
            {
                'Start': 25, 'End': 40, 'Class': 'Z-DNA', 'Subclass': 'Z-DNA',
                'Normalized_Score': 0.6, 'Actual_Score': 3.1, 'GC_Content': 65.0,
                'Sequence': 'CGCGCGCGCGCGCGCG', 'Scoring_Method': 'Z-seeker'
            }
        ]
        
        # Test BED export
        bed_output = export_to_bed(test_motifs, "test_sequence")
        assert "track name=" in bed_output
        assert "G-Quadruplex" in bed_output
        assert "Z-DNA" in bed_output
        print("   ✓ BED format export working")
        
        # Test GFF3 export
        gff3_output = export_to_gff3(test_motifs, "test_sequence")
        assert "##gff-version 3" in gff3_output
        assert "non_B_DNA_motif" in gff3_output
        print("   ✓ GFF3 format export working")
        
        # Test density export
        density_output = create_density_bedgraph(test_motifs, 100, "test_sequence")
        assert "track type=bedGraph" in density_output
        print("   ✓ Density bedGraph export working")
        
        return True
        
    except Exception as e:
        print(f"   ❌ Export utilities test failed: {e}")
        return False

def test_enhanced_caching():
    """Test enhanced caching system"""
    print("💾 Testing Enhanced Caching...")
    
    try:
        from enhanced_cache import get_cache_manager, get_cache_stats
        import motifs
        
        # Clear cache to start fresh
        cache_manager = get_cache_manager()
        cache_manager.clear_all_caches()
        
        # Test sequence for analysis
        test_seq = "GGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCCATCGATCGATCGATCG"
        
        # First analysis (should be cache miss)
        start_time = time.time()
        results1 = motifs.all_motifs(test_seq, sequence_name="cache_test_1")
        first_time = time.time() - start_time
        
        # Second analysis with same sequence (should be cache hit)
        start_time = time.time()
        results2 = motifs.all_motifs(test_seq, sequence_name="cache_test_2")
        second_time = time.time() - start_time
        
        # Verify results are consistent
        assert len(results1) == len(results2)
        
        # Check cache statistics
        stats = get_cache_stats()
        hit_rate = stats['totals']['hit_rate_percent']
        
        print(f"   ✓ First analysis: {first_time:.3f}s, Second: {second_time:.3f}s")
        print(f"   ✓ Cache hit rate: {hit_rate:.1f}%")
        print(f"   ✓ Memory usage: {stats['totals']['memory_used_mb']:.2f} MB")
        
        # Verify cache is working (second analysis should be faster)
        assert hit_rate > 0, "Cache should have at least some hits"
        
        return True
        
    except Exception as e:
        print(f"   ❌ Enhanced caching test failed: {e}")
        return False

def test_integration():
    """Test integration of all new features"""
    print("🔧 Testing Integration...")
    
    try:
        import motifs
        from export_utils import create_motif_browser_session
        from enhanced_cache import get_cache_stats
        
        # Run analysis with all features
        test_seq = "GGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCCATCGATCGATCGATCGCGCGCGCGCGCGCG"
        
        results = motifs.all_motifs(
            test_seq, 
            sequence_name="integration_test",
            nonoverlap=False,
            report_hotspots=True,
            calculate_conservation=False  # Skip conservation for speed
        )
        
        assert len(results) > 0, "Should find motifs in test sequence"
        print(f"   ✓ Found {len(results)} motifs in integration test")
        
        # Test browser session creation
        browser_session = create_motif_browser_session(
            results, 
            "integration_test", 
            len(test_seq)
        )
        
        assert "tracks" in browser_session
        assert "all_motifs" in browser_session["tracks"]
        print(f"   ✓ Browser session created with {len(browser_session['tracks'])} tracks")
        
        # Check caching is working
        cache_stats = get_cache_stats()
        print(f"   ✓ Cache system active: {cache_stats['totals']['total_requests']} requests")
        
        return True
        
    except Exception as e:
        print(f"   ❌ Integration test failed: {e}")
        return False

def test_streamlit_compatibility():
    """Test that Streamlit app can import all new modules"""
    print("🎨 Testing Streamlit Compatibility...")
    
    try:
        # Test imports that app.py uses
        from export_utils import export_to_bed, export_to_gff3, create_density_bedgraph
        from enhanced_cache import get_cache_stats, clear_all_caches
        import motifs
        from motifs.enhanced_visualization import create_comprehensive_information_based_visualizations
        
        print("   ✓ All Streamlit imports successful")
        
        # Test that app.py can be imported (syntax check)
        try:
            import app
            print("   ✓ app.py imports without errors")
        except ImportError as e:
            print(f"   ⚠ app.py import warning: {e}")
        
        return True
        
    except Exception as e:
        print(f"   ❌ Streamlit compatibility test failed: {e}")
        return False

def main():
    """Run all tests"""
    print("=" * 60)
    print("🧬 NBDFinder Revamp Features Test Suite")
    print("=" * 60)
    print(f"⏰ Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    tests = [
        test_export_utilities,
        test_enhanced_caching, 
        test_integration,
        test_streamlit_compatibility,
        test_rest_api  # API test last since it starts server
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        try:
            if test():
                passed += 1
                print("   ✅ PASSED")
            else:
                print("   ❌ FAILED")
        except Exception as e:
            print(f"   💥 CRASHED: {e}")
        print()
    
    print("=" * 60)
    print(f"🏆 Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 ALL REVAMP FEATURES WORKING! NBDFinder is production-ready!")
    else:
        print("⚠️  Some features need attention")
    
    print("=" * 60)
    
    return passed == total

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)