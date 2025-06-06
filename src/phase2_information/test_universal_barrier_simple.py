"""
Simple unit tests for universal barrier - no external dependencies.

This can run immediately to verify our Step 1 code audit is complete.
"""

def test_pattern_analyzer_barrier():
    """Test that PatternAnalyzer always returns 2 coins."""
    # Import here to avoid module-level dependency
    from pattern_analyzer import PatternAnalyzer
    
    analyzer = PatternAnalyzer()
    
    # Test various inputs - barrier should always be 2
    test_cases = [
        (set(), 1),
        (set([0]), 1),
        (set([0, 1, 2, 3, 4, 5, 6, 7]), 1),
        (set([2, 3]), 5),
        (set(range(100)), 100),
    ]
    
    print("Testing PatternAnalyzer._calculate_barrier()...")
    for unique_rungs, n_components in test_cases:
        barrier = analyzer._calculate_barrier(unique_rungs, n_components)
        assert barrier == 2, f"Expected 2 coins, got {barrier}"
        print(f"  ✓ rungs={len(unique_rungs)}, components={n_components} → {barrier} coins")
    
    print("PatternAnalyzer: ✅ Always returns 2 coins\n")

def test_pattern_analyzer_v2_barrier():
    """Test that PatternAnalyzerV2 always returns 0.18 eV."""
    from pattern_analyzer_v2 import PatternAnalyzerV2, E_COH
    
    analyzer = PatternAnalyzerV2()
    
    # Test various phase properties - barrier should always be 0.18 eV
    test_cases = [
        (0.1, 0.001, 1.0, 0.0),
        (100.0, 1000.0, 10000.0, 1.0),
        (5.0, 10.0, 100.0, 0.5),
        (20.0, 50.0, 500.0, 0.8),
    ]
    
    print("Testing PatternAnalyzerV2._calculate_effective_barrier()...")
    expected = 2 * E_COH
    for coherence, flow, density, frustration in test_cases:
        barrier = analyzer._calculate_effective_barrier(
            coherence_length=coherence,
            info_flow=flow,
            rec_density=density,
            frustration=frustration
        )
        assert abs(barrier - expected) < 1e-10, f"Expected {expected} eV, got {barrier}"
        print(f"  ✓ coherence={coherence}, flow={flow} → {barrier:.3f} eV")
    
    print(f"PatternAnalyzerV2: ✅ Always returns {expected:.3f} eV\n")

def verify_constants():
    """Verify that E_COH is consistent across modules."""
    from pattern_analyzer import E_COH as E_COH_V1
    from pattern_analyzer_v2 import E_COH as E_COH_V2
    
    print("Verifying E_COH consistency...")
    print(f"  PatternAnalyzer: E_COH = {E_COH_V1} eV")
    print(f"  PatternAnalyzerV2: E_COH = {E_COH_V2} eV")
    
    assert E_COH_V1 == E_COH_V2, "E_COH mismatch between modules!"
    assert E_COH_V1 == 0.090, f"E_COH = {E_COH_V1}, expected 0.090 eV"
    
    print("  ✅ E_COH = 0.090 eV in both modules\n")

def main():
    """Run all simple tests."""
    print("=" * 60)
    print("Universal Barrier Simple Tests")
    print("=" * 60)
    print()
    
    # Verify constants
    verify_constants()
    
    # Test barrier calculations
    test_pattern_analyzer_barrier()
    test_pattern_analyzer_v2_barrier()
    
    print("=" * 60)
    print("✅ ALL TESTS PASSED!")
    print("Recognition Science principles confirmed:")
    print("  - Barrier is always 2 coins (0.18 eV)")
    print("  - No protein-specific variations")
    print("  - E_COH = 0.090 eV universally")
    print("=" * 60)

if __name__ == "__main__":
    main() 