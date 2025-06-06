#!/usr/bin/env python3
"""Debug full protein template analysis."""

import numpy as np
from accelerated_folder_v3 import AcceleratedFolderV3

# Test with Trp-cage
sequence = 'NLYIQWLKDGGPSSGRPPPS'
n_residues = len(sequence)

print("Debugging Trp-cage Template Analysis")
print("="*60)

# Create folder
folder = AcceleratedFolderV3(
    n_residues=n_residues,
    sequence=sequence,
    temperature=296,
    use_template_params=True
)

# Run template formation
print("\nPhase 1: Template formation...")
template_metrics = folder._run_until_template_with_analysis()

if template_metrics['template_formed']:
    print(f"Template formed at {template_metrics['template_time_us']:.3f} μs")
    
    # Get the analysis
    analysis = folder.template_analysis
    
    print(f"\nVoxel Graph Analysis:")
    print(f"  N voxels: {analysis.n_voxels}")
    print(f"  N components: {analysis.n_components}")
    print(f"  N loops: {analysis.n_loops}")
    print(f"  Contact order: {analysis.contact_order:.3f}")
    
    print(f"\nSecondary Structure:")
    print(f"  Helix fraction: {analysis.helix_fraction:.3f}")
    print(f"  Sheet fraction: {analysis.sheet_fraction:.3f}")
    
    print(f"\nRS-Derived Metrics:")
    print(f"  Path entropy: {analysis.path_entropy:.6f}")
    print(f"  Mobility anisotropy: {analysis.mobility_anisotropy:.6f}")
    
    print(f"\nCalculated Parameters:")
    print(f"  Barrier: {analysis.barrier_coins} coins ({analysis.barrier_ev} eV)")
    print(f"  P_ledger: {analysis.p_ledger:.6f}")
    print(f"  P_geom: {analysis.p_geom:.6f}")
    
    # Check if P_geom is φ^(-1)
    PHI = 1.618033988749895
    if abs(analysis.p_geom - 1/PHI) < 1e-6:
        print(f"  ⚠️  P_geom = φ^(-1) exactly!")
    
    # Calculate what P_geom should be based on the analysis
    expected_p_geom = PHI ** (-analysis.n_loops / 2)
    if analysis.contact_order > 0:
        expected_p_geom *= np.exp(-analysis.contact_order)
    if analysis.helix_fraction > 0.3:
        expected_p_geom *= PHI ** (analysis.helix_fraction / 2)
    if analysis.sheet_fraction > 0.3:
        expected_p_geom *= PHI ** (-analysis.sheet_fraction)
    expected_p_geom = np.clip(expected_p_geom, 1e-4, 1.0)
    
    print(f"\n  Expected P_geom from components: {expected_p_geom:.6f}")
    
    # Calculate final k0 modifier
    path_efficiency = 1.0 - analysis.path_entropy
    mobility_factor = 1.0 / (1.0 + analysis.mobility_anisotropy)
    k0_modifier = analysis.p_ledger * analysis.p_geom * path_efficiency * mobility_factor
    
    print(f"\nFinal k0 calculation:")
    print(f"  Path efficiency (1 - entropy): {path_efficiency:.6f}")
    print(f"  Mobility factor: {mobility_factor:.6f}")
    print(f"  Total k0 modifier: {k0_modifier:.6f}")
    
    # Compare to default
    default_modifier = 0.5 * 0.01  # p_ledger=0.5, p_geom=0.01
    print(f"  Default modifier: {default_modifier:.6f}")
    print(f"  Ratio: {k0_modifier/default_modifier:.2f}x")
else:
    print("Template did not form!") 