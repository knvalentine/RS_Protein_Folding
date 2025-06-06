"""
Physical-Layer k₀ Calibration from RS First Principles

This module implements the correct k₀(N, topology) formula based on
Recognition Science axioms. The barrier is always 0.18 eV (2 coins),
and only the prefactor varies with protein properties.
"""

# RS Fundamental Constants
E_COH = 0.090  # eV (one recognition quantum)
TAU_0 = 7.33e-15  # s (fundamental tick)
PHI = 1.618033988749895  # Golden ratio

def calculate_k0_rs_pure(n_residues: int, p_ledger: float = 0.5, 
                        p_geom: float = 0.01) -> float:
    """
    Calculate k₀ from Recognition Science first principles.
    
    k₀ = (1/8τ₀) × Φ^(-N/2) × P_ledger × P_geom
    
    This is the Arrhenius prefactor. The barrier is always 2 coins (0.18 eV).
    
    Args:
        n_residues: Number of residues (N)
        p_ledger: Ledger availability factor (0-1)
        p_geom: Geometric compatibility factor (0-1)
        
    Returns:
        k₀ in s^-1
    """
    # Base rate: eight-beat cycle frequency
    base_rate = 1 / (8 * TAU_0)
    
    # Size-dependent phase alignment probability
    # Larger proteins need more residues to be phase-coherent
    phase_factor = PHI ** (-n_residues / 2)
    
    # Total prefactor
    k0 = base_rate * phase_factor * p_ledger * p_geom
    
    return k0

def calculate_folding_rate(k0: float, temperature_k: float = 310.0) -> float:
    """
    Calculate actual folding rate using Arrhenius equation.
    
    k_fold = k₀ × exp(-ΔC/kT)
    
    where ΔC = 0.18 eV (universal RS barrier)
    
    Args:
        k0: Prefactor in s^-1
        temperature_k: Temperature in Kelvin
        
    Returns:
        k_fold in s^-1
    """
    kT = 8.617e-5 * temperature_k  # eV
    barrier = 0.18  # eV (2 coins, universal)
    
    k_fold = k0 * np.exp(-barrier / kT)
    return k_fold

def calculate_mean_folding_time(k0: float, temperature_k: float = 310.0) -> float:
    """
    Calculate mean folding time.
    
    <t> = 1/k_fold = 1/(k₀ × exp(-ΔC/kT))
    
    Args:
        k0: Prefactor in s^-1
        temperature_k: Temperature in Kelvin
        
    Returns:
        Mean folding time in seconds
    """
    k_fold = calculate_folding_rate(k0, temperature_k)
    return 1 / k_fold

def show_k0_scaling():
    """
    Demonstrate how k₀ scales with protein size.
    """
    print("=" * 60)
    print("k₀ Scaling with Protein Size (RS First Principles)")
    print("=" * 60)
    print("\nFormula: k₀ = (1/8τ₀) × Φ^(-N/2) × P_ledger × P_geom")
    print(f"Base rate: 1/8τ₀ = {1/(8*TAU_0):.2e} s^-1")
    print(f"P_ledger = 0.5 (default)")
    print(f"P_geom = 0.01 (default)")
    print(f"Barrier = 0.18 eV (universal)\n")
    
    print("N    k₀ (s^-1)      <t> at 310K")
    print("-" * 40)
    
    for n in [5, 10, 20, 30, 40, 50]:
        k0 = calculate_k0_rs_pure(n)
        mean_time = calculate_mean_folding_time(k0, 310.0)
        mean_time_us = mean_time * 1e6
        
        print(f"{n:2d}   {k0:8.2e}   {mean_time_us:10.2f} μs")
    
    print("\nNote: These are default values. Real proteins will have")
    print("different P_ledger and P_geom based on their topology.")

def calibrate_from_template(template_analysis):
    """
    Calculate k₀ using parameters extracted from template.
    
    Args:
        template_analysis: TemplateAnalysis object from pattern_analyzer
        
    Returns:
        Dictionary with k₀ and related parameters
    """
    # Extract parameters from template
    n_residues = template_analysis.n_voxels  # Approximate
    p_ledger = template_analysis.p_ledger
    p_geom = template_analysis.p_geom
    
    # Calculate k₀
    k0 = calculate_k0_rs_pure(n_residues, p_ledger, p_geom)
    
    # Calculate rates at different temperatures
    temps = [280, 300, 310, 320, 340]
    rates = {}
    times = {}
    
    for T in temps:
        k_fold = calculate_folding_rate(k0, T)
        rates[T] = k_fold
        times[T] = 1 / k_fold
    
    return {
        'k0': k0,
        'n_residues': n_residues,
        'p_ledger': p_ledger,
        'p_geom': p_geom,
        'barrier_ev': 0.18,  # Always!
        'rates_by_temp': rates,
        'times_by_temp': times
    }

def validate_against_experiments():
    """
    Show how our k₀ values compare to experimental data.
    """
    print("\n" + "=" * 60)
    print("Validation Against Experimental Data")
    print("=" * 60)
    
    # Known experimental folding times
    experiments = [
        ("Trp-cage", 20, 4.1e-6, 296),  # 4.1 μs at 296K
        ("WW domain", 34, 13.0e-6, 298),  # 13 μs at 298K
        ("Villin", 35, 0.7e-6, 300),  # 0.7 μs at 300K
        ("BBA5", 23, 13.0e-6, 298),  # 13 μs at 298K
    ]
    
    print("\nProtein     N   Exp. time   Temp   RS k₀      RS time    Ratio")
    print("-" * 70)
    
    for name, n, exp_time, temp in experiments:
        # Use default parameters for now
        k0 = calculate_k0_rs_pure(n, p_ledger=0.5, p_geom=0.01)
        rs_time = calculate_mean_folding_time(k0, temp)
        ratio = rs_time / exp_time
        
        print(f"{name:10s} {n:2d}  {exp_time*1e6:6.1f} μs  {temp}K  "
              f"{k0:8.2e}  {rs_time*1e6:6.1f} μs  {ratio:6.2f}")
    
    print("\nNote: Using default P_ledger=0.5, P_geom=0.01")
    print("Template-specific values will improve agreement.")

# Import numpy only if available
try:
    import numpy as np
except ImportError:
    # Define exp function without numpy
    import math
    class np:
        @staticmethod
        def exp(x):
            return math.exp(x)

if __name__ == "__main__":
    # Show k₀ scaling
    show_k0_scaling()
    
    # Validate against experiments
    validate_against_experiments()
    
    print("\n" + "=" * 60)
    print("Key Points:")
    print("- Barrier is ALWAYS 0.18 eV (2 coins)")
    print("- Only k₀ varies with protein size/topology")
    print("- Template analysis provides P_ledger and P_geom")
    print("- No hidden scaling factors!")
    print("=" * 60) 