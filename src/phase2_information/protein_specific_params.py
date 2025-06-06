"""
Protein-Specific Parameter Calculator for Recognition Science

This module integrates:
1. Torsion analysis → barrier height
2. Geometric factors → P_geom
3. Ledger factors → P_ledger (simplified for now)

All parameters derived from first principles based on
Deeper Understanding.txt guidelines.
"""

import numpy as np
from typing import Dict, Tuple

from torsion_analysis import get_torsion_summary
from geometric_factor import calculate_geometric_factor

# Golden ratio
PHI = (1 + np.sqrt(5)) / 2

# RS Constants
TAU_0 = 7.33e-15  # s
E_COH = 0.090  # eV


def calculate_ledger_factor(length: int, n_domains: int) -> float:
    """
    Calculate ledger availability factor P_ledger.
    
    Based on Section 7 of Deeper Understanding:
    - Ledger occupancy f ≈ 3.3 × 10^-22
    - P_ledger = φ^(-N_free/2) where N_free is degrees of freedom
    
    Args:
        length: Number of residues
        n_domains: Number of structural domains
        
    Returns:
        P_ledger factor
    """
    # Degrees of freedom for ledger balancing
    # Fewer domains = more freedom = higher P_ledger
    n_free = max(1, n_domains - 1)
    
    # Base ledger factor
    p_ledger = PHI ** (-n_free / 2)
    
    # Size correction (smaller proteins have higher ledger availability)
    size_factor = min(1.5, 30 / length)
    
    return p_ledger * size_factor


def get_protein_specific_params(sequence: str) -> Dict:
    """
    Calculate all protein-specific parameters from sequence.
    
    Args:
        sequence: Amino acid sequence
        
    Returns:
        Dictionary with all parameters needed for folding simulation
    """
    # Get torsion analysis
    torsion = get_torsion_summary(sequence)
    
    # Get geometric factor
    geom = calculate_geometric_factor(
        sequence,
        torsion['helix_content'],
        torsion['sheet_content'],
        torsion['domains']
    )
    
    # Get ledger factor
    p_ledger = calculate_ledger_factor(
        torsion['length'],
        len(torsion['domains'])
    )
    
    # Calculate size-dependent base prefactor
    # k₀(n) = (1/8τ₀) × φ^(-n/2)
    base_k0 = (1 / (8 * TAU_0)) * PHI ** (-torsion['length'] / 2)
    
    # Total prefactor
    k0_total = base_k0 * p_ledger * geom['p_geom']
    
    # Compile results
    params = {
        'sequence': sequence,
        'length': torsion['length'],
        
        # Barrier
        'barrier_coins': torsion['barrier_coins'],
        'barrier_ev': torsion['barrier_ev'],
        
        # Prefactor components
        'base_k0': base_k0,
        'p_ledger': p_ledger,
        'p_geom': geom['p_geom'],
        'k0_total': k0_total,
        
        # Structure info
        'structure_type': geom['structure_type'],
        'helix_content': torsion['helix_content'],
        'sheet_content': torsion['sheet_content'],
        'n_domains': len(torsion['domains']),
        
        # Detailed components for analysis
        'torsion_summary': torsion,
        'geometric_details': geom
    }
    
    return params


def predict_folding_time(params: Dict, temperature: float = 298.0) -> Dict:
    """
    Predict folding time using protein-specific parameters.
    
    Args:
        params: Output from get_protein_specific_params
        temperature: Temperature in K
        
    Returns:
        Folding time predictions
    """
    kT = 8.617e-5 * temperature  # eV
    
    # Folding rate
    k_fold = params['k0_total'] * np.exp(-params['barrier_ev'] / kT)
    
    # Mean folding time
    mean_time_s = 1 / k_fold
    mean_time_us = mean_time_s * 1e6
    
    return {
        'k_fold': k_fold,
        'mean_time_s': mean_time_s,
        'mean_time_us': mean_time_us,
        'temperature_k': temperature
    }


# Test and compare with our fixed parameters
if __name__ == "__main__":
    test_proteins = {
        'Trp-cage': {
            'sequence': 'NLYIQWLKDGGPSSGRPPPS',
            'exp_us': 4.1,
            'exp_temp': 296
        },
        'Villin': {
            'sequence': 'LSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF',
            'exp_us': 0.7,
            'exp_temp': 300
        },
        'BBA5': {
            'sequence': 'EQYTAKYKGRTFRNEKELRDFIE',
            'exp_us': 13.0,
            'exp_temp': 298
        },
        'WW domain': {
            'sequence': 'GSKLPPGWEKRMSRSSGRVYYFNHITNASQWERPS',
            'exp_us': 13.0,
            'exp_temp': 298
        }
    }
    
    print("Protein-Specific Parameter Analysis")
    print("=" * 70)
    
    for name, protein in test_proteins.items():
        params = get_protein_specific_params(protein['sequence'])
        prediction = predict_folding_time(params, protein['exp_temp'])
        
        print(f"\n{name}:")
        print(f"  Barrier: {params['barrier_coins']} coins ({params['barrier_ev']:.3f} eV)")
        print(f"  vs fixed: 2 coins (0.180 eV)")
        print(f"  Prefactor components:")
        print(f"    Base k₀: {params['base_k0']:.2e} s⁻¹")
        print(f"    P_ledger: {params['p_ledger']:.3f}")
        print(f"    P_geom: {params['p_geom']:.4f}")
        print(f"    Total k₀: {params['k0_total']:.2e} s⁻¹")
        
        # Compare with fixed parameters
        fixed_k0 = params['base_k0'] * 0.5 * 0.01  # Our fixed values
        print(f"  vs fixed k₀: {fixed_k0:.2e} s⁻¹")
        print(f"  Ratio: {params['k0_total']/fixed_k0:.1f}x")
        
        print(f"\n  Folding time prediction:")
        print(f"    Protein-specific: {prediction['mean_time_us']:.1f} μs")
        print(f"    Experimental: {protein['exp_us']:.1f} μs")
        print(f"    Ratio (RS/Exp): {prediction['mean_time_us']/protein['exp_us']:.2f}")
        
        # What would fixed params predict?
        k_fixed = fixed_k0 * np.exp(-0.18 / (8.617e-5 * protein['exp_temp']))
        fixed_time_us = 1e6 / k_fixed
        print(f"    Fixed params would give: {fixed_time_us:.1f} μs")
        print(f"    Fixed ratio (RS/Exp): {fixed_time_us/protein['exp_us']:.2f}") 