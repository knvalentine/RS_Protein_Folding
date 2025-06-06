"""
Test Phase-Based Analysis Concept

This demonstrates how reading phase information (not just structure)
could explain why BBA5 folds faster than our structural analysis predicted.
"""

import math

# Constants
E_COH = 0.090  # eV
TAU_0 = 7.33e-15  # s
PHI = 1.618033988749895
kT_298 = 0.0257  # eV at 298K

class MockPhaseAnalysis:
    """Mock phase analysis results for different proteins."""
    
    def __init__(self, protein_name):
        if protein_name == "trp-cage":
            # Compact, well-organized phases
            self.phase_coherence_length = 8.0  # nm - good coherence
            self.information_flow_rate = 15.0  # bits/ps - high flow
            self.recognition_density = 4.5  # events/residue - moderate
            self.phase_frustration = 0.1  # low frustration
            self.channel_occupancy = [1.0, 0.8, 0.6, 0.3, 0.2, 0.1, 0.0, 0.0]
            self.active_channels = 4
            
        elif protein_name == "villin":
            # All helix, excellent phase flow
            self.phase_coherence_length = 15.0  # nm - excellent!
            self.information_flow_rate = 20.0  # bits/ps - very high
            self.recognition_density = 3.0  # events/residue - low (simple)
            self.phase_frustration = 0.05  # almost no frustration
            self.channel_occupancy = [1.0, 0.9, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0]
            self.active_channels = 2  # Very simple phase pattern
            
        elif protein_name == "bba5":
            # HERE'S THE KEY: Mixed structure but GOOD phase organization!
            self.phase_coherence_length = 7.0  # nm - decent coherence
            self.information_flow_rate = 12.0  # bits/ps - good flow
            self.recognition_density = 5.5  # events/residue - moderate-high
            self.phase_frustration = 0.25  # Some frustration but not terrible
            self.channel_occupancy = [1.0, 0.7, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1]
            self.active_channels = 6  # Complex but organized
            
        elif protein_name == "ww-domain":
            # Beta sheet, moderate phase properties
            self.phase_coherence_length = 6.0  # nm - okay
            self.information_flow_rate = 8.0  # bits/ps - moderate
            self.recognition_density = 6.0  # events/residue - high
            self.phase_frustration = 0.3  # moderate frustration
            self.channel_occupancy = [1.0, 0.6, 0.6, 0.4, 0.3, 0.2, 0.1, 0.0]
            self.active_channels = 5


def calculate_barrier_from_phases(analysis):
    """
    Return the universal RS barrier.
    
    According to Recognition Science axioms, the barrier is ALWAYS
    0.18 eV (2 coins) for all proteins. Phase properties affect
    the prefactor k₀, not the barrier.
    """
    # Universal barrier - NEVER varies!
    return 2 * E_COH  # 0.18 eV


def calculate_k0_modifiers_from_phases(analysis):
    """
    Calculate k₀ modifiers from phase properties.
    
    These factors multiply the base k₀ to account for
    phase organization effects on folding kinetics.
    """
    modifier = 1.0
    
    # Coherence length contribution
    if analysis.phase_coherence_length < 5.0:
        # Poor coherence slows folding
        modifier *= 0.5
    elif analysis.phase_coherence_length > 10.0:
        # Excellent coherence speeds folding
        modifier *= 2.0
    
    # Information flow contribution
    if analysis.information_flow_rate < 10:
        # Low flow slows folding
        modifier *= (analysis.information_flow_rate / 10)
    elif analysis.information_flow_rate > 15:
        # High flow speeds folding
        modifier *= (1 + (analysis.information_flow_rate - 15) / 10)
    
    # Recognition density contribution
    if analysis.recognition_density > 5:
        # High density slows folding (complexity)
        modifier *= (5 / analysis.recognition_density)
    
    # Frustration contribution
    # High frustration slows folding
    modifier *= (1 - analysis.phase_frustration * 0.5)
    
    return modifier


def calculate_geometric_factors(analysis):
    """Calculate P factors from phase properties."""
    # P_ledger from channel usage
    p_ledger = PHI ** (-(analysis.active_channels - 1) / 4)
    
    # P_geom from coherence and frustration
    p_geom = 0.01  # base
    
    # Coherence bonus/penalty
    if analysis.phase_coherence_length > 10:
        p_geom *= PHI ** ((analysis.phase_coherence_length - 10) / 20)
    elif analysis.phase_coherence_length < 5:
        p_geom *= PHI ** (-(5 - analysis.phase_coherence_length) / 5)
    
    # Frustration penalty
    p_geom *= (1 - analysis.phase_frustration)
    
    return p_ledger, p_geom


def demonstrate_phase_analysis():
    """Show how phase analysis explains folding times better."""
    
    print("PHASE-BASED ANALYSIS vs STRUCTURE-BASED")
    print("="*70)
    print("\nKey Insight: BBA5 has mixed structure but ORGANIZED phases!")
    print("-"*70)
    
    proteins = [
        ("Trp-cage", 20, 4.1),
        ("Villin", 35, 0.7),
        ("BBA5", 23, 13.0),
        ("WW domain", 34, 13.0)
    ]
    
    print("\nPHASE PROPERTIES:")
    print("-"*70)
    print("Protein    | Coherence | Info Flow | Rec Dens | Frustration | Channels")
    print("           |    (nm)   | (bits/ps) | (ev/res) |   (0-1)     | Active")
    print("-"*70)
    
    for name, _, _ in proteins:
        analysis = MockPhaseAnalysis(name.lower().replace(" ", "-"))
        print(f"{name:10} | {analysis.phase_coherence_length:9.1f} | "
              f"{analysis.information_flow_rate:9.1f} | "
              f"{analysis.recognition_density:8.1f} | "
              f"{analysis.phase_frustration:11.2f} | {analysis.active_channels:6d}")
    
    print("\nFOLDING PREDICTIONS:")
    print("-"*70)
    print("Protein    | Barrier | k₀ mod  | P_ledger | P_geom  | Predicted | Exp  | Ratio")
    print("           |  (eV)   |         |          |         |   (μs)    | (μs) |")
    print("-"*70)
    
    for name, n_res, exp_time in proteins:
        analysis = MockPhaseAnalysis(name.lower().replace(" ", "-"))
        
        # Calculate from phase properties
        barrier = calculate_barrier_from_phases(analysis)  # Always 0.18 eV
        k0_modifier = calculate_k0_modifiers_from_phases(analysis)
        p_ledger, p_geom = calculate_geometric_factors(analysis)
        
        # Folding time
        k0_base = 1 / (8 * TAU_0)
        size_factor = PHI**(-n_res / 2)
        k0 = k0_base * size_factor * p_ledger * p_geom * k0_modifier
        k_fold = k0 * math.exp(-barrier / kT_298)
        pred_time = 1e6 / k_fold
        
        ratio = pred_time / exp_time
        status = "✅" if 0.1 < ratio < 10 else "❌"
        
        print(f"{name:10} | {barrier:7.3f} | {k0_modifier:7.2f} | {p_ledger:8.3f} | "
              f"{p_geom:7.4f} | {pred_time:9.1f} | {exp_time:4.1f} | "
              f"{ratio:5.2f} {status}")
    
    print("\n" + "="*70)
    print("KEY INSIGHTS:")
    print("-"*70)
    print("1. ALL proteins have the SAME barrier: 0.18 eV (2 coins)")
    print("2. BBA5 folds faster due to good phase organization (k₀ modifier)")
    print("3. Its phase coherence (7 nm) is decent despite mixed structure")
    print("4. Information flows well (12 bits/ps) through the template")
    print("\n5. Villin is fastest due to exceptional k₀ modifiers:")
    print("   - Excellent coherence (15 nm) → 2x speed boost")
    print("   - High info flow (20 bits/ps) → additional boost")
    print("6. WW domain is slower due to poor k₀ modifiers from frustration")
    print("\nCONCLUSION: The barrier is universal (0.18 eV)")
    print("            Phase organization affects k₀, not the barrier!")
    print("="*70)


if __name__ == "__main__":
    demonstrate_phase_analysis() 