#!/usr/bin/env python3
"""Check what parameters the template analysis extracts."""

from accelerated_folder_v3 import AcceleratedFolderV3
import numpy as np

# Test proteins
test_cases = [
    ('Trp-cage', 'NLYIQWLKDGGPSSGRPPPS', 296),
    ('BBA5', 'EQYTAKYKGRTFRNEKELRDFIE', 298),
    ('Villin', 'LSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF', 300)
]

print("Template-Extracted Parameters Analysis")
print("="*60)

for name, sequence, temp in test_cases:
    print(f"\n{name} ({len(sequence)} residues)")
    print("-"*40)
    
    folder = AcceleratedFolderV3(
        n_residues=len(sequence),
        temperature=temp,
        sequence=sequence,
        use_template_params=True
    )
    
    # Run to get parameters
    result = folder.run_accelerated()
    
    if 'parameters_used' in result:
        params = result['parameters_used']
        print(f"Barrier: {params['barrier_ev']:.3f} eV (coins: {params['barrier_coins']})")
        print(f"P_ledger: {params['p_ledger']:.3f}")
        print(f"P_geom: {params['p_geom']:.4f}")
        print(f"Path entropy: {params['path_entropy']:.4f}")
        print(f"Mobility anisotropy: {params['mobility_anisotropy']:.4f}")
        print(f"N_voxels: {params['n_voxels']}")
        print(f"N_components: {params['n_components']}")
        print(f"N_loops: {params['n_loops']}")
        print(f"Contact order: {params['contact_order']:.3f}")
        
        # Calculate effective k0 modifier
        k0_modifier = params['p_ledger'] * params['p_geom'] * params['path_entropy'] * (1/(1+params['mobility_anisotropy']))
        print(f"\nEffective k0 modifier: {k0_modifier:.6f}")
        print(f"  (compared to default: p_ledger=0.5, p_geom=0.01 → 0.005)")
        print(f"  Ratio to default: {k0_modifier/0.005:.2f}x") 