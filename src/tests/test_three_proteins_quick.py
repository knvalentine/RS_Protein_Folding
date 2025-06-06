#!/usr/bin/env python3
"""Quick test of 3 proteins to verify setup."""

import numpy as np
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'core'))
from accelerated_folder_v3 import AcceleratedFolderV3

# Test just 3 proteins
test_proteins = [
    ('Trp-cage', 'NLYIQWLKDGGPSSGRPPPS', 4.1, 296),
    ('BBA5', 'EQYTAKYKGRTFRNEKELRDFIE', 13.0, 298),
    ('Villin', 'LSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF', 0.7, 300)
]

print("Quick 3-Protein Test")
print("="*60)

np.random.seed(42)

for name, sequence, exp_time, temp in test_proteins:
    print(f"\n{name} ({len(sequence)} residues)")
    print("-"*40)
    
    folder = AcceleratedFolderV3(
        n_residues=len(sequence),
        sequence=sequence,
        temperature=temp,
        use_template_params=True
    )
    
    # Single run
    result = folder.run_accelerated(max_us=10000.0)
    
    if result.get('total_time_us'):
        total_time = result['total_time_us']
        ratio = total_time / exp_time
        
        print(f"RS prediction: {total_time:.1f} μs")
        print(f"Experimental: {exp_time} μs")
        print(f"Ratio: {ratio:.2f}")
        
        if 0.5 <= ratio <= 2.0:
            print("✅ Within factor of 2!")
        elif 0.1 <= ratio <= 10.0:
            print("⚠️  Within order of magnitude")
        else:
            print(f"❌ Off by {max(ratio, 1/ratio):.0f}x")

print("\nQuick test complete!") 