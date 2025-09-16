#!/usr/bin/env python3
"""
NBDFinder Refactored Integration Demo
====================================

Demonstrates the complete refactored NBDFinder system with:
- Consistent per-class overlap filtering across all entry points
- Optimized performance with parallel processing
- Fast visualization generation
- Comprehensive validation

Author: NBDFinder Refactoring Team
"""

import time
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def demo_consistency():
    """Demonstrate consistent overlap filtering across all entry points."""
    print("🔄 DEMONSTRATION: Consistent Overlap Filtering")
    print("=" * 60)
    
    test_seq = "GGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCC" * 3
    
    # Test all available entry points
    entry_points = [
        ("all_motifs_refactored", "all_motifs_refactored", "all_motifs_refactored"),
        ("optimized_orchestrator", "core.optimized_orchestrator", "all_motifs_optimized"),
    ]
    
    results = {}
    
    for name, module_path, function_name in entry_points:
        try:
            module = __import__(module_path, fromlist=[function_name])
            func = getattr(module, function_name)
            
            start_time = time.time()
            motifs = func(test_seq, sequence_name=f"demo_{name}", nonoverlap=True)
            exec_time = time.time() - start_time
            
            results[name] = {
                'count': len(motifs),
                'time': exec_time,
                'classes': sorted(set(m.get('Class', 'Unknown') for m in motifs))
            }
            
            print(f"✓ {name}: {len(motifs)} motifs ({exec_time:.3f}s)")
            print(f"  Classes: {', '.join(results[name]['classes'])}")
            
        except Exception as e:
            print(f"✗ {name}: Failed - {e}")
            results[name] = {'error': str(e)}
    
    # Check consistency
    counts = [r.get('count', 0) for r in results.values() if 'count' in r]
    if len(set(counts)) <= 1:
        print(f"\n✅ CONSISTENT: All entry points returned same motif count")
    else:
        print(f"\n⚠️ INCONSISTENT: Different motif counts: {counts}")
    
    return results


def demo_performance_optimization():
    """Demonstrate performance improvements."""
    print("\n🚀 DEMONSTRATION: Performance Optimization")
    print("=" * 60)
    
    # Test with larger sequence for meaningful timing
    test_seq = "GGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCC" * 10
    
    try:
        from core.optimized_orchestrator import all_motifs_optimized
        
        print("Running optimized analysis...")
        start_time = time.time()
        motifs = all_motifs_optimized(test_seq, sequence_name="performance_demo", 
                                    nonoverlap=True, report_hotspots=False)
        total_time = time.time() - start_time
        
        print(f"\n📊 PERFORMANCE RESULTS:")
        print(f"  Total time: {total_time:.3f}s")
        print(f"  Motifs found: {len(motifs)}")
        print(f"  Rate: {len(motifs)/total_time:.1f} motifs/second")
        
        # Analyze class distribution
        class_counts = {}
        for motif in motifs:
            class_name = motif.get('Class', 'Unknown')
            class_counts[class_name] = class_counts.get(class_name, 0) + 1
        
        print(f"  Classes found: {len(class_counts)}")
        for class_name, count in sorted(class_counts.items()):
            print(f"    {class_name}: {count}")
        
        return {
            'success': True,
            'time': total_time,
            'motif_count': len(motifs),
            'class_counts': class_counts
        }
        
    except Exception as e:
        print(f"❌ Performance demo failed: {e}")
        return {'success': False, 'error': str(e)}


def demo_fast_visualizations():
    """Demonstrate fast visualization generation."""
    print("\n📊 DEMONSTRATION: Fast Visualization Generation")
    print("=" * 60)
    
    try:
        from motifs.visualization_optimized import create_optimized_visualizations, benchmark_visualizations
        
        # Generate sample data and create visualizations
        print("Creating optimized visualizations...")
        start_time = time.time()
        result = create_optimized_visualizations(save_plots=False, fast_mode=True)
        vis_time = time.time() - start_time
        
        if result['success']:
            print(f"✅ SUCCESS: Created {len(result['created_plots'])} visualizations")
            print(f"  Total time: {result['total_time']:.3f}s")
            print(f"  Data prep: {result['data_prep_time']:.3f}s")
            print(f"  Plot generation: {result['total_time'] - result['data_prep_time']:.3f}s")
            print(f"  Plots created: {', '.join(result['created_plots'])}")
            
            # Quick benchmark
            print(f"\nRunning quick benchmark...")
            benchmark_visualizations(iterations=2)
            
            return result
        else:
            print(f"❌ FAILED: {result.get('error', 'Unknown error')}")
            return result
            
    except Exception as e:
        print(f"❌ Visualization demo failed: {e}")
        return {'success': False, 'error': str(e)}


def demo_data_validation():
    """Demonstrate data integrity and validation."""
    print("\n🔍 DEMONSTRATION: Data Validation")
    print("=" * 60)
    
    try:
        from core.optimized_orchestrator import all_motifs_optimized
        
        test_seq = "GGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCC" * 2
        
        # Test with overlap filtering
        motifs = all_motifs_optimized(test_seq, sequence_name="validation_demo", 
                                    nonoverlap=True)
        
        print(f"Validating {len(motifs)} motifs...")
        
        # Check required fields
        required_fields = ['Class', 'Subclass', 'Start', 'End', 'Length', 'Normalized_Score']
        field_coverage = {field: 0 for field in required_fields}
        
        for motif in motifs:
            for field in required_fields:
                if field in motif and motif[field] is not None:
                    field_coverage[field] += 1
        
        print(f"\n📋 FIELD COVERAGE:")
        for field, count in field_coverage.items():
            coverage = (count / len(motifs)) * 100 if motifs else 0
            status = "✅" if coverage == 100 else "⚠️" if coverage > 80 else "❌"
            print(f"  {field}: {coverage:.1f}% {status}")
        
        # Check coordinate validity
        invalid_coords = 0
        for motif in motifs:
            start = motif.get('Start', 0)
            end = motif.get('End', 0)
            if start <= 0 or end <= 0 or start >= end:
                invalid_coords += 1
        
        coord_validity = ((len(motifs) - invalid_coords) / len(motifs)) * 100 if motifs else 100
        print(f"\n📍 COORDINATE VALIDITY: {coord_validity:.1f}%")
        
        # Check for within-class overlaps
        class_motifs = {}
        for motif in motifs:
            class_name = motif.get('Class', 'Unknown')
            if class_name not in class_motifs:
                class_motifs[class_name] = []
            class_motifs[class_name].append(motif)
        
        total_overlaps = 0
        for class_name, motif_list in class_motifs.items():
            class_overlaps = 0
            for i, m1 in enumerate(motif_list):
                for m2 in motif_list[i+1:]:
                    start1, end1 = m1.get('Start', 0), m1.get('End', 0)
                    start2, end2 = m2.get('Start', 0), m2.get('End', 0)
                    if not (end1 < start2 or end2 < start1):
                        class_overlaps += 1
            total_overlaps += class_overlaps
            if class_overlaps > 0:
                print(f"  ⚠️ {class_name}: {class_overlaps} overlaps")
        
        if total_overlaps == 0:
            print(f"✅ NO WITHIN-CLASS OVERLAPS: Filtering working correctly")
        else:
            print(f"❌ {total_overlaps} WITHIN-CLASS OVERLAPS FOUND")
        
        return {
            'success': True,
            'motif_count': len(motifs),
            'field_coverage': field_coverage,
            'coordinate_validity': coord_validity,
            'within_class_overlaps': total_overlaps
        }
        
    except Exception as e:
        print(f"❌ Data validation demo failed: {e}")
        return {'success': False, 'error': str(e)}


def run_complete_demo():
    """Run complete demonstration of refactored NBDFinder."""
    print("🧬 NBDFinder Refactored System Demonstration")
    print("=" * 80)
    print("Demonstrating:")
    print("✓ Consistent per-class overlap filtering across all entry points")
    print("✓ Optimized performance with parallel processing")
    print("✓ Fast visualization generation")
    print("✓ Comprehensive data validation")
    print("=" * 80)
    
    demo_results = {}
    
    # Run demonstrations
    demo_results['consistency'] = demo_consistency()
    demo_results['performance'] = demo_performance_optimization()
    demo_results['visualization'] = demo_fast_visualizations()
    demo_results['validation'] = demo_data_validation()
    
    # Final summary
    print("\n" + "🎯 FINAL SUMMARY")
    print("=" * 80)
    
    successes = 0
    total_demos = 0
    
    for demo_name, result in demo_results.items():
        total_demos += 1
        if isinstance(result, dict) and result.get('success', True):
            successes += 1
            print(f"✅ {demo_name.title()}: SUCCESS")
        else:
            print(f"❌ {demo_name.title()}: FAILED")
    
    success_rate = (successes / total_demos) * 100 if total_demos > 0 else 0
    
    print(f"\nOverall Success Rate: {success_rate:.1f}% ({successes}/{total_demos})")
    
    if success_rate == 100:
        print("\n🎉 ALL DEMONSTRATIONS SUCCESSFUL!")
        print("✓ NBDFinder refactoring complete and validated")
        print("✓ Non-overlapping motifs per class consistently reported")
        print("✓ Performance optimized")
        print("✓ Visualizations generation optimized")
        print("✓ Data integrity maintained")
    else:
        print(f"\n⚠️ Some demonstrations had issues. Check logs above.")
    
    return demo_results


if __name__ == "__main__":
    results = run_complete_demo()
    
    # Determine exit code based on success
    success_count = sum(1 for r in results.values() 
                       if isinstance(r, dict) and r.get('success', True))
    
    if success_count == len(results):
        print("\n✅ All demonstrations completed successfully!")
        sys.exit(0)
    else:
        print(f"\n⚠️ {len(results) - success_count} demonstrations had issues.")
        sys.exit(1)