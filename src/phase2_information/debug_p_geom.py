#!/usr/bin/env python3
"""Debug P_geom calculation."""

import numpy as np
from pattern_analyzer import PatternAnalyzer

PHI = 1.618033988749895

# Test the P_geom calculation directly
analyzer = PatternAnalyzer()

# Test cases
test_cases = [
    {"n_loops": 0, "contact_order": 0.0, "helix_frac": 0.0, "sheet_frac": 0.0},
    {"n_loops": 0, "contact_order": 0.0, "helix_frac": 0.5, "sheet_frac": 0.0},
    {"n_loops": 0, "contact_order": 0.0, "helix_frac": 0.0, "sheet_frac": 0.5},
    {"n_loops": 1, "contact_order": 0.1, "helix_frac": 0.3, "sheet_frac": 0.0},
    {"n_loops": 2, "contact_order": 0.2, "helix_frac": 0.0, "sheet_frac": 0.3},
]

print("P_geom Calculation Debug")
print("="*60)

for i, case in enumerate(test_cases):
    print(f"\nTest case {i+1}:")
    print(f"  n_loops: {case['n_loops']}")
    print(f"  contact_order: {case['contact_order']}")
    print(f"  helix_frac: {case['helix_frac']}")
    print(f"  sheet_frac: {case['sheet_frac']}")
    
    # Calculate step by step
    p_geom = PHI ** (-case['n_loops'] / 2)
    print(f"  After loops: {p_geom:.6f}")
    
    if case['contact_order'] > 0:
        p_geom *= np.exp(-case['contact_order'])
        print(f"  After contact order: {p_geom:.6f}")
    
    if case['helix_frac'] > 0.3:
        p_geom *= PHI ** (case['helix_frac'] / 2)
        print(f"  After helix bonus: {p_geom:.6f}")
    
    if case['sheet_frac'] > 0.3:
        p_geom *= PHI ** (-case['sheet_frac'])
        print(f"  After sheet penalty: {p_geom:.6f}")
    
    # Final clipping
    p_geom_final = np.clip(p_geom, 1e-4, 1.0)
    print(f"  Final (after clip): {p_geom_final:.6f}")
    
    # Check if it's φ^(-1)
    if abs(p_geom_final - 1/PHI) < 1e-6:
        print("  ⚠️  This equals φ^(-1) = 0.618034!")

# Check what φ^(-1) actually is
print(f"\n\nφ^(-1) = {1/PHI:.6f}")
print(f"φ^(-1/2) = {PHI**(-0.5):.6f}") 