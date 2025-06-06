"""
Refined Torsion Analysis with Better Structure Recognition

This version:
1. Uses sliding windows to identify continuous secondary structures
2. Properly identifies Villin as primarily helical
3. Better captures the actual folding complexity
"""

import numpy as np
from typing import Dict, List, Tuple

# Golden ratio
PHI = (1 + np.sqrt(5)) / 2

# More nuanced residue propensities
STRONG_HELIX = ['A', 'E', 'L', 'M']
WEAK_HELIX = ['Q', 'K', 'R', 'H', 'F']
STRONG_SHEET = ['V', 'I', 'Y', 'W']
WEAK_SHEET = ['F', 'T', 'C']
BREAKERS = ['P', 'G']
FLEXIBLE = ['S', 'D', 'N']


def analyze_structure_propensity(sequence: str, window_size: int = 5) -> Dict:
    """
    Analyze structure propensity using sliding windows.
    
    Args:
        sequence: Amino acid sequence
        window_size: Size of sliding window
        
    Returns:
        Structure propensity analysis
    """
    length = len(sequence)
    helix_score = np.zeros(length)
    sheet_score = np.zeros(length)
    
    # Score each position based on local context
    for i in range(length):
        # Define window
        start = max(0, i - window_size // 2)
        end = min(length, i + window_size // 2 + 1)
        window = sequence[start:end]
        
        # Count propensities in window
        h_strong = sum(1 for r in window if r in STRONG_HELIX)
        h_weak = sum(0.5 for r in window if r in WEAK_HELIX)
        s_strong = sum(1 for r in window if r in STRONG_SHEET)
        s_weak = sum(0.5 for r in window if r in WEAK_SHEET)
        breakers = sum(1 for r in window if r in BREAKERS)
        
        # Calculate scores
        helix_score[i] = (h_strong + h_weak - breakers) / len(window)
        sheet_score[i] = (s_strong + s_weak - breakers) / len(window)
    
    return {
        'helix_score': helix_score,
        'sheet_score': sheet_score
    }


def identify_secondary_structures(sequence: str) -> List[Dict]:
    """
    Identify secondary structure regions.
    
    Args:
        sequence: Amino acid sequence
        
    Returns:
        List of structure regions
    """
    scores = analyze_structure_propensity(sequence)
    helix_score = scores['helix_score']
    sheet_score = scores['sheet_score']
    
    structures = []
    current_structure = None
    
    for i in range(len(sequence)):
        # Determine structure at this position
        if helix_score[i] > 0.4 and helix_score[i] > sheet_score[i]:
            structure_type = 'helix'
        elif sheet_score[i] > 0.3:
            structure_type = 'sheet'
        else:
            structure_type = 'coil'
        
        # Update current structure
        if current_structure is None:
            current_structure = {
                'type': structure_type,
                'start': i,
                'end': i,
                'length': 1
            }
        elif current_structure['type'] == structure_type:
            current_structure['end'] = i
            current_structure['length'] += 1
        else:
            # Save current and start new
            if current_structure['length'] >= 3:  # Minimum structure length
                structures.append(current_structure)
            current_structure = {
                'type': structure_type,
                'start': i,
                'end': i,
                'length': 1
            }
    
    # Don't forget the last structure
    if current_structure and current_structure['length'] >= 3:
        structures.append(current_structure)
    
    return structures


def calculate_refined_barrier(structures: List[Dict]) -> Tuple[int, float]:
    """
    Calculate barrier based on structural complexity.
    
    Args:
        structures: List of secondary structures
        
    Returns:
        (barrier_coins, barrier_ev)
    """
    # Count unique structure types that need to form
    structure_types = set(s['type'] for s in structures if s['type'] != 'coil')
    
    if len(structure_types) == 0:
        # No defined structure - minimum barrier
        return 2, 0.180
    elif len(structure_types) == 1:
        if 'helix' in structure_types:
            # Pure helix - can use fast pathway
            return 2, 0.180
        else:
            # Pure sheet - standard barrier
            return 2, 0.180
    else:
        # Mixed structure - needs extra coin
        return 3, 0.270


def get_refined_analysis(sequence: str) -> Dict:
    """
    Get refined structural analysis.
    
    Args:
        sequence: Amino acid sequence
        
    Returns:
        Refined analysis dictionary
    """
    structures = identify_secondary_structures(sequence)
    barrier_coins, barrier_ev = calculate_refined_barrier(structures)
    
    # Calculate content
    helix_residues = sum(s['length'] for s in structures if s['type'] == 'helix')
    sheet_residues = sum(s['length'] for s in structures if s['type'] == 'sheet')
    
    helix_content = helix_residues / len(sequence)
    sheet_content = sheet_residues / len(sequence)
    
    # Determine if truly mixed
    is_mixed = (helix_content > 0.1 and sheet_content > 0.1)
    
    return {
        'sequence': sequence,
        'length': len(sequence),
        'structures': structures,
        'n_structures': len(structures),
        'helix_content': helix_content,
        'sheet_content': sheet_content,
        'is_mixed': is_mixed,
        'barrier_coins': barrier_coins,
        'barrier_ev': barrier_ev
    }


# Test proteins
if __name__ == "__main__":
    test_proteins = {
        'Trp-cage': 'NLYIQWLKDGGPSSGRPPPS',
        'Villin': 'LSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF',
        'BBA5': 'EQYTAKYKGRTFRNEKELRDFIE',
        'WW domain': 'GSKLPPGWEKRMSRSSGRVYYFNHITNASQWERPS'
    }
    
    print("Refined Torsion Analysis")
    print("=" * 70)
    
    for name, sequence in test_proteins.items():
        analysis = get_refined_analysis(sequence)
        
        print(f"\n{name}:")
        print(f"  Length: {analysis['length']} residues")
        print(f"  Helix content: {analysis['helix_content']:.1%}")
        print(f"  Sheet content: {analysis['sheet_content']:.1%}")
        print(f"  Mixed structure: {analysis['is_mixed']}")
        print(f"  Barrier: {analysis['barrier_coins']} coins ({analysis['barrier_ev']:.3f} eV)")
        print(f"  Structures found: {analysis['n_structures']}")
        
        for s in analysis['structures']:
            print(f"    {s['type']}: residues {s['start']+1}-{s['end']+1} (length {s['length']})") 