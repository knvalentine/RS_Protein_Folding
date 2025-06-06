"""
Torsion Analysis for Recognition Science Protein Folding

This module analyzes protein sequences to:
1. Assign expected φ/ψ torsion rungs for each residue
2. Count coin costs based on deviation from golden ratio geometry
3. Identify simultaneous rung locks needed

Based on Deeper Understanding.txt:
- α-helix: 100° ≈ φ² × 60° base → +1 coin
- β-sheet: 180° ≈ φ³ × 60° base → +2 coins
- Golden ratio minimizes J(θ) = ½(θ + 1/θ)
"""

import numpy as np
from typing import Dict, List, Tuple, Optional

# Golden ratio
PHI = (1 + np.sqrt(5)) / 2

# Torsion rung assignments based on secondary structure propensities
# From statistical analysis of PDB structures
HELIX_FORMERS = ['A', 'E', 'L', 'M', 'Q', 'K', 'R', 'H']
SHEET_FORMERS = ['V', 'I', 'Y', 'W', 'F', 'T', 'C']
TURN_FORMERS = ['G', 'P', 'S', 'D', 'N']

# Torsion states and their coin costs
TORSION_STATES = {
    'helix': {
        'phi': -60,  # degrees
        'psi': -45,
        'coin_cost': 1,  # φ² deviation
        'rung': 32  # Base rung for helix
    },
    'sheet': {
        'phi': -120,
        'psi': 120,
        'coin_cost': 2,  # φ³ deviation
        'rung': 33  # Higher rung for sheet
    },
    'turn': {
        'phi': -60,
        'psi': -30,
        'coin_cost': 0,  # Golden ratio aligned
        'rung': 31  # Lower rung for turns
    },
    'coil': {
        'phi': -60,
        'psi': -40,
        'coin_cost': 0,  # Flexible, no lock needed
        'rung': 30  # Base rung
    }
}


def analyze_sequence(sequence: str) -> Dict:
    """
    Analyze protein sequence to assign torsion states and coin costs.
    
    Args:
        sequence: Amino acid sequence (single letter code)
        
    Returns:
        Dictionary with torsion analysis results
    """
    residue_states = []
    total_coin_cost = 0
    rung_distribution = {}
    
    # Assign torsion state to each residue based on propensity
    for i, residue in enumerate(sequence):
        if residue in HELIX_FORMERS:
            state = 'helix'
        elif residue in SHEET_FORMERS:
            state = 'sheet'
        elif residue in TURN_FORMERS:
            state = 'turn'
        else:
            state = 'coil'
        
        state_info = TORSION_STATES[state]
        residue_states.append({
            'residue': residue,
            'position': i,
            'state': state,
            'phi': state_info['phi'],
            'psi': state_info['psi'],
            'coin_cost': state_info['coin_cost'],
            'rung': state_info['rung']
        })
        
        total_coin_cost += state_info['coin_cost']
        rung = state_info['rung']
        rung_distribution[rung] = rung_distribution.get(rung, 0) + 1
    
    # Identify domains and simultaneous locks needed
    domains = identify_domains(residue_states)
    simultaneous_locks = count_simultaneous_locks(domains)
    
    return {
        'sequence': sequence,
        'length': len(sequence),
        'residue_states': residue_states,
        'total_coin_cost': total_coin_cost,
        'rung_distribution': rung_distribution,
        'domains': domains,
        'simultaneous_locks': simultaneous_locks,
        'helix_content': sum(1 for r in residue_states if r['state'] == 'helix') / len(sequence),
        'sheet_content': sum(1 for r in residue_states if r['state'] == 'sheet') / len(sequence),
        'mixed_structure': ('helix' in [r['state'] for r in residue_states] and 
                           'sheet' in [r['state'] for r in residue_states])
    }


def identify_domains(residue_states: List[Dict]) -> List[Dict]:
    """
    Identify continuous domains of secondary structure.
    
    Args:
        residue_states: List of residue state dictionaries
        
    Returns:
        List of domain dictionaries
    """
    domains = []
    current_domain = None
    
    for i, res_state in enumerate(residue_states):
        state = res_state['state']
        
        # Skip coil regions
        if state == 'coil':
            if current_domain:
                domains.append(current_domain)
                current_domain = None
            continue
        
        # Start new domain or extend current
        if not current_domain or current_domain['type'] != state:
            if current_domain:
                domains.append(current_domain)
            current_domain = {
                'type': state,
                'start': i,
                'end': i,
                'length': 1,
                'rung': res_state['rung'],
                'coin_cost': res_state['coin_cost']
            }
        else:
            current_domain['end'] = i
            current_domain['length'] += 1
    
    if current_domain:
        domains.append(current_domain)
    
    return domains


def count_simultaneous_locks(domains: List[Dict]) -> int:
    """
    Count number of distinct rungs that must lock simultaneously.
    
    This determines the barrier height:
    - Single rung type: 2 coins (minimum barrier)
    - Multiple rungs: 2 + extra coins for each additional rung
    
    Args:
        domains: List of domain dictionaries
        
    Returns:
        Number of simultaneous locks needed
    """
    unique_rungs = set()
    for domain in domains:
        # Only count structured domains that need locking
        if domain['type'] in ['helix', 'sheet']:
            unique_rungs.add(domain['rung'])
    
    return len(unique_rungs)


def calculate_barrier_coins(analysis: Dict) -> int:
    """
    Calculate total barrier height in coins.
    
    Based on Deeper Understanding.txt:
    - Minimum 2 coins (leave unfolded + enter folded)
    - +1 coin for each additional simultaneous rung lock
    
    Args:
        analysis: Results from analyze_sequence
        
    Returns:
        Barrier height in number of coins
    """
    base_coins = 2  # Minimum for any folding transition
    extra_coins = max(0, analysis['simultaneous_locks'] - 1)
    
    return base_coins + extra_coins


def get_torsion_summary(sequence: str) -> Dict:
    """
    Get complete torsion analysis summary for a protein sequence.
    
    Args:
        sequence: Amino acid sequence
        
    Returns:
        Summary dictionary with all relevant parameters
    """
    analysis = analyze_sequence(sequence)
    barrier_coins = calculate_barrier_coins(analysis)
    barrier_ev = barrier_coins * 0.090  # E_coh = 0.090 eV
    
    summary = {
        'sequence': sequence,
        'length': analysis['length'],
        'helix_content': analysis['helix_content'],
        'sheet_content': analysis['sheet_content'],
        'mixed_structure': analysis['mixed_structure'],
        'total_coin_cost': analysis['total_coin_cost'],
        'barrier_coins': barrier_coins,
        'barrier_ev': barrier_ev,
        'simultaneous_locks': analysis['simultaneous_locks'],
        'domains': analysis['domains'],
        'rung_distribution': analysis['rung_distribution']
    }
    
    return summary


# Test with our proteins
if __name__ == "__main__":
    test_proteins = {
        'Trp-cage': 'NLYIQWLKDGGPSSGRPPPS',
        'Villin': 'LSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF',
        'BBA5': 'EQYTAKYKGRTFRNEKELRDFIE',
        'WW domain': 'GSKLPPGWEKRMSRSSGRVYYFNHITNASQWERPS'
    }
    
    print("Torsion Analysis Results")
    print("=" * 70)
    
    for name, sequence in test_proteins.items():
        summary = get_torsion_summary(sequence)
        print(f"\n{name}:")
        print(f"  Length: {summary['length']} residues")
        print(f"  Helix content: {summary['helix_content']:.1%}")
        print(f"  Sheet content: {summary['sheet_content']:.1%}")
        print(f"  Mixed structure: {summary['mixed_structure']}")
        print(f"  Barrier: {summary['barrier_coins']} coins ({summary['barrier_ev']:.3f} eV)")
        print(f"  Simultaneous locks: {summary['simultaneous_locks']}")
        print(f"  Domains: {len(summary['domains'])}")
        for domain in summary['domains']:
            print(f"    - {domain['type']}: residues {domain['start']+1}-{domain['end']+1}") 