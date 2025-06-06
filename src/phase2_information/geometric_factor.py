"""
Geometric Factor Calculation for Recognition Science

This module calculates P_geom based on:
1. Contact order of the protein
2. Voxel walk topology
3. Loop-sum contributions (Section 13 of Deeper Understanding)

P_geom represents the fraction of voxel paths that are geometrically
compatible with the folded structure.
"""

import numpy as np
from typing import Dict, List, Tuple

# Golden ratio
PHI = (1 + np.sqrt(5)) / 2

# Voxel parameters
VOXEL_SIZE = 0.335  # nm (from Deeper Understanding.txt)


def calculate_contact_order(sequence: str, structure_type: str) -> float:
    """
    Calculate relative contact order based on structure type.
    
    Contact order = average sequence separation of contacting residues / length
    
    Args:
        sequence: Amino acid sequence
        structure_type: Dominant structure ('helix', 'sheet', 'mixed')
        
    Returns:
        Relative contact order (0-1)
    """
    length = len(sequence)
    
    if structure_type == 'helix':
        # Helices have mostly local contacts (i, i+3 and i, i+4)
        avg_separation = 3.6  # Average for helix
        return avg_separation / length
    
    elif structure_type == 'sheet':
        # Sheets can have long-range contacts
        avg_separation = length / 3  # Rough estimate
        return min(1.0, avg_separation / length)
    
    else:  # mixed
        # Intermediate between helix and sheet
        avg_separation = length / 5
        return min(1.0, avg_separation / length)


def calculate_loop_degeneracy(length: int, contact_order: float) -> float:
    """
    Calculate loop degeneracy factor from voxel walk statistics.
    
    Based on Section 13: Loop-sum Σ₁, Σ₂, ... give partition function
    
    Args:
        length: Number of residues
        contact_order: Relative contact order
        
    Returns:
        Loop degeneracy factor
    """
    # Number of voxel transitions
    n_transitions = length / PHI  # From Deeper Understanding
    
    # Effective loop length
    effective_loop = length * contact_order
    
    # Loop degeneracy scales as φ^(-loop_length/2)
    # This comes from the damping factor A = √P · φ^(-1/2) for each edge
    loop_factor = PHI ** (-effective_loop / 2)
    
    return loop_factor


def calculate_helix_bonus(helix_content: float) -> float:
    """
    Calculate geometric bonus for helical proteins.
    
    Helices have optimal golden-ratio geometry, so they get a bonus.
    
    Args:
        helix_content: Fraction of residues in helix
        
    Returns:
        Helix bonus factor (>= 1)
    """
    # Helices align with golden ratio → bonus factor
    # Maximum bonus of φ for pure helix
    return 1 + (PHI - 1) * helix_content


def calculate_sheet_penalty(sheet_content: float) -> float:
    """
    Calculate geometric penalty for sheet-rich proteins.
    
    Sheets require more complex topology → penalty.
    
    Args:
        sheet_content: Fraction of residues in sheet
        
    Returns:
        Sheet penalty factor (<= 1)
    """
    # Sheets have φ³ deviation → harder to form
    # Penalty scales with sheet content
    return 1 - 0.5 * sheet_content


def get_structure_type(helix_content: float, sheet_content: float) -> str:
    """
    Determine dominant structure type.
    
    Args:
        helix_content: Fraction helix
        sheet_content: Fraction sheet
        
    Returns:
        'helix', 'sheet', or 'mixed'
    """
    if helix_content > 0.6 and sheet_content < 0.2:
        return 'helix'
    elif sheet_content > 0.4:
        return 'sheet'
    else:
        return 'mixed'


def calculate_geometric_factor(sequence: str, helix_content: float, 
                             sheet_content: float, domains: List[Dict]) -> Dict:
    """
    Calculate complete geometric factor P_geom.
    
    Args:
        sequence: Amino acid sequence
        helix_content: Fraction of helix
        sheet_content: Fraction of sheet
        domains: List of domain dictionaries from torsion analysis
        
    Returns:
        Dictionary with P_geom and components
    """
    length = len(sequence)
    structure_type = get_structure_type(helix_content, sheet_content)
    
    # Base factors
    contact_order = calculate_contact_order(sequence, structure_type)
    loop_degeneracy = calculate_loop_degeneracy(length, contact_order)
    
    # Structure-specific modifiers
    helix_bonus = calculate_helix_bonus(helix_content)
    sheet_penalty = calculate_sheet_penalty(sheet_content)
    
    # Domain complexity penalty
    # More domains = more complex folding topology
    domain_penalty = PHI ** (-len(domains) / 10)
    
    # Combine all factors
    p_geom = loop_degeneracy * helix_bonus * sheet_penalty * domain_penalty
    
    # Ensure probability-like bounds [1e-8, 1.0]
    p_geom = max(1e-8, min(1.0, p_geom))
    
    return {
        'p_geom': p_geom,
        'structure_type': structure_type,
        'contact_order': contact_order,
        'loop_degeneracy': loop_degeneracy,
        'helix_bonus': helix_bonus,
        'sheet_penalty': sheet_penalty,
        'domain_penalty': domain_penalty,
        'n_domains': len(domains)
    }


# Test with our proteins
if __name__ == "__main__":
    from torsion_analysis import get_torsion_summary
    
    test_proteins = {
        'Trp-cage': 'NLYIQWLKDGGPSSGRPPPS',
        'Villin': 'LSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF',
        'BBA5': 'EQYTAKYKGRTFRNEKELRDFIE',
        'WW domain': 'GSKLPPGWEKRMSRSSGRVYYFNHITNASQWERPS'
    }
    
    print("Geometric Factor Analysis")
    print("=" * 70)
    
    for name, sequence in test_proteins.items():
        # Get torsion analysis first
        torsion = get_torsion_summary(sequence)
        
        # Calculate geometric factor
        geom = calculate_geometric_factor(
            sequence,
            torsion['helix_content'],
            torsion['sheet_content'],
            torsion['domains']
        )
        
        print(f"\n{name}:")
        print(f"  Structure type: {geom['structure_type']}")
        print(f"  P_geom: {geom['p_geom']:.4f}")
        print(f"  Components:")
        print(f"    Contact order: {geom['contact_order']:.3f}")
        print(f"    Loop degeneracy: {geom['loop_degeneracy']:.4f}")
        print(f"    Helix bonus: {geom['helix_bonus']:.3f}")
        print(f"    Sheet penalty: {geom['sheet_penalty']:.3f}")
        print(f"    Domain penalty: {geom['domain_penalty']:.3f}")
        
        # Compare to our fixed value
        print(f"  vs fixed P_geom: {0.01:.4f}")
        print(f"  Ratio: {geom['p_geom']/0.01:.1f}x") 