#!/usr/bin/env python3
"""
Demonstration of Hyperscan Registry Loader Integration

This script demonstrates:
1. Registry loading from precompiled files
2. Automatic fallback when Hyperscan is unavailable
3. Scanner initialization with registry preloading
4. Detection with merged 10-mer regions
"""

import os
import sys

# Ensure we can import from the repo
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def demo_registry_loading():
    """Demonstrate registry loading capabilities"""
    print("=" * 70)
    print("DEMO 1: Registry Loading")
    print("=" * 70)
    
    from utils.motif_patterns import get_hs_db_for_class, get_pattern_registry
    
    # Load Z-DNA registry
    print("\n[Loading Z-DNA registry]")
    db, id_to_ten, id_to_score = get_hs_db_for_class('ZDNA', 'registry')
    
    print(f"  • Database object: {type(db).__name__ if db else 'None (fallback mode)'}")
    print(f"  • Number of patterns: {len(id_to_ten)}")
    print(f"  • Score range: {min(id_to_score.values()):.2f} - {max(id_to_score.values()):.2f}")
    print(f"  • Sample patterns: {list(id_to_ten.values())[:3]}")
    
    # Load A-philic registry
    print("\n[Loading A-philic registry]")
    db, id_to_ten, id_to_score = get_hs_db_for_class('APhilic', 'registry')
    
    print(f"  • Database object: {type(db).__name__ if db else 'None (fallback mode)'}")
    print(f"  • Number of patterns: {len(id_to_ten)}")
    print(f"  • Score range: {min(id_to_score.values()):.2f} - {max(id_to_score.values()):.2f}")
    print(f"  • Sample patterns: {list(id_to_ten.values())[:3]}")
    
    # Show registry metadata
    print("\n[Registry Metadata]")
    registry = get_pattern_registry('ZDNA', 'registry')
    print(f"  • Class: {registry['class']}")
    print(f"  • Generated: {registry['generated_at']}")
    print(f"  • Source: {registry['meta']['source']}")


def demo_scanner_integration():
    """Demonstrate scanner integration"""
    print("\n" + "=" * 70)
    print("DEMO 2: Scanner Integration")
    print("=" * 70)
    
    from utils.modular_scanner import ModularMotifDetector
    
    # Initialize scanner (automatically preloads registries)
    print("\n[Initializing Scanner]")
    scanner = ModularMotifDetector(registry_dir='registry')
    
    print(f"  ✓ Scanner initialized")
    if hasattr(scanner, 'hsdb_map'):
        print(f"  ✓ Preloaded registries for: {', '.join(scanner.hsdb_map.keys())}")
        for cls_name, info in scanner.hsdb_map.items():
            n_patterns = len(info['id_to_ten'])
            print(f"    - {cls_name}: {n_patterns} patterns")


def demo_detection():
    """Demonstrate motif detection with registry-based matching"""
    print("\n" + "=" * 70)
    print("DEMO 3: Motif Detection")
    print("=" * 70)
    
    from utils.modular_scanner import ModularMotifDetector
    
    # Create test sequences
    test_sequences = {
        'Z-DNA': 'ATCGCGCGCGCGCGCGCGCGAT',  # CG repeats
        'A-philic': 'AGGGGGGGGGAGGGGGGGGC' + 'AAAA' + 'CCCCCCCCCCCCCCCCCCCC',  # G/C-rich
        'Mixed': 'CGCGCGCG' + 'AAAAA' + 'AGGGGGGGGG' + 'TTTTT' + 'GCGCGCGC',
    }
    
    scanner = ModularMotifDetector(registry_dir='registry')
    
    for name, seq in test_sequences.items():
        print(f"\n[Testing {name} sequence]")
        print(f"  Sequence: {seq[:50]}{'...' if len(seq) > 50 else ''}")
        print(f"  Length: {len(seq)} bp")
        
        motifs = scanner.analyze_sequence(seq, name)
        
        if motifs:
            print(f"  ✓ Found {len(motifs)} motifs:")
            for motif in motifs:
                cls = motif.get('Class', 'Unknown')
                subcls = motif.get('Subclass', '')
                start = motif.get('Start', 0)
                end = motif.get('End', 0)
                score = motif.get('Score', 0)
                print(f"    - {cls}/{subcls}: [{start:3d}-{end:3d}] Score: {score:.3f}")
        else:
            print(f"  No motifs detected")


def demo_detector_attributes():
    """Demonstrate detector class attributes"""
    print("\n" + "=" * 70)
    print("DEMO 4: Detector Class Attributes")
    print("=" * 70)
    
    from motif_detection.z_dna_detector import ZDNADetector
    from motif_detection.a_philic_detector import APhilicDetector
    from utils.modular_scanner import ModularMotifDetector
    
    # Initialize scanner to set class attributes
    print("\n[Initializing Scanner to Set Attributes]")
    scanner = ModularMotifDetector(registry_dir='registry')
    
    # Check Z-DNA detector
    print("\n[ZDNADetector]")
    if hasattr(ZDNADetector, 'HS_DB_INFO'):
        info = ZDNADetector.HS_DB_INFO
        print(f"  ✓ HS_DB_INFO is set")
        print(f"  • Contains {len(info['id_to_ten'])} patterns")
        print(f"  • Database: {type(info['db']).__name__ if info['db'] else 'None'}")
    else:
        print(f"  HS_DB_INFO not set (optional)")
    
    # Check A-philic detector
    print("\n[APhilicDetector]")
    if hasattr(APhilicDetector, 'HS_DB_INFO'):
        info = APhilicDetector.HS_DB_INFO
        print(f"  ✓ HS_DB_INFO is set")
        print(f"  • Contains {len(info['id_to_ten'])} patterns")
        print(f"  • Database: {type(info['db']).__name__ if info['db'] else 'None'}")
    else:
        print(f"  HS_DB_INFO not set (optional)")


def demo_caching():
    """Demonstrate in-memory caching"""
    print("\n" + "=" * 70)
    print("DEMO 5: In-Memory Caching")
    print("=" * 70)
    
    from utils.motif_patterns import get_hs_db_for_class
    import time
    
    print("\n[First Load - Cache Miss]")
    start = time.time()
    db1, id_to_ten1, id_to_score1 = get_hs_db_for_class('ZDNA', 'registry')
    time1 = (time.time() - start) * 1000
    print(f"  Time: {time1:.2f} ms")
    
    print("\n[Second Load - Cache Hit]")
    start = time.time()
    db2, id_to_ten2, id_to_score2 = get_hs_db_for_class('ZDNA', 'registry')
    time2 = (time.time() - start) * 1000
    print(f"  Time: {time2:.2f} ms")
    
    print(f"\n  ✓ Same object returned: {id_to_ten1 is id_to_ten2}")
    print(f"  • Speedup: {time1/time2:.1f}x faster")


def main():
    """Run all demonstrations"""
    print("\n" + "=" * 70)
    print("HYPERSCAN REGISTRY LOADER INTEGRATION DEMO")
    print("=" * 70)
    
    try:
        demo_registry_loading()
        demo_scanner_integration()
        demo_detection()
        demo_detector_attributes()
        demo_caching()
        
        print("\n" + "=" * 70)
        print("DEMO COMPLETE ✓")
        print("=" * 70)
        
        print("\nKey Features Demonstrated:")
        print("  • Registry loading from precompiled files")
        print("  • Automatic Hyperscan fallback when unavailable")
        print("  • Scanner initialization with registry preloading")
        print("  • Motif detection with merged 10-mer regions")
        print("  • Detector class attribute population")
        print("  • In-memory caching for performance")
        
    except Exception as e:
        print(f"\n[ERROR] Demo failed: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
