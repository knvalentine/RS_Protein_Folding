"""
Final test of protein-specific parameters showing actual improvements
"""

import numpy as np
from accelerated_folder_v2 import AcceleratedFolderV2

def test_protein_final(name: str, sequence: str, exp_time_us: float, exp_temp_k: float):
    """Test with clear display of parameters used."""
    
    print(f"\n{'='*70}")
    print(f"{name} - Final Test")
    print(f"{'='*70}")
    print(f"Experimental: {exp_time_us} μs at {exp_temp_k}K")
    
    # Default parameters
    folder_default = AcceleratedFolderV2(
        n_residues=len(sequence),
        temperature=exp_temp_k,
        sequence=sequence,
        use_protein_specific=False
    )
    
    # Get one sample to show parameters
    time_us, _ = folder_default.monte_carlo_barrier_crossing()
    
    print(f"\nDEFAULT PARAMETERS:")
    print(f"  Barrier: 0.180 eV (2 coins)")
    print(f"  k₀: {folder_default.k0_used:.2e} s⁻¹")
    print(f"  k(T): {folder_default.k_fold:.2e} s⁻¹")
    print(f"  Sample time: {time_us:.1f} μs")
    print(f"  Ratio: {time_us/exp_time_us:.2f}")
    
    # Protein-specific parameters
    folder_specific = AcceleratedFolderV2(
        n_residues=len(sequence),
        temperature=exp_temp_k,
        sequence=sequence,
        use_protein_specific=True
    )
    
    params = folder_specific.get_parameter_summary()
    
    # Get one sample
    time_us_specific, _ = folder_specific.monte_carlo_barrier_crossing()
    
    print(f"\nPROTEIN-SPECIFIC PARAMETERS:")
    print(f"  Structure: {params['helix_content']:.0%} helix, {params['sheet_content']:.0%} sheet")
    print(f"  Barrier: {params['barrier_ev']:.3f} eV ({params['barrier_coins']} coins)")
    print(f"  P_ledger: {params['p_ledger']:.2f}")
    print(f"  P_geom: {params['p_geom']:.3f}")
    print(f"  k₀: {folder_specific.k0_used:.2e} s⁻¹")
    print(f"  k(T): {folder_specific.k_fold:.2e} s⁻¹")
    print(f"  Sample time: {time_us_specific:.1f} μs")
    print(f"  Ratio: {time_us_specific/exp_time_us:.2f}")
    
    # Run multiple samples for statistics
    print(f"\nSTATISTICS (10 samples each):")
    
    default_times = []
    for _ in range(10):
        t, _ = folder_default.monte_carlo_barrier_crossing()
        default_times.append(t)
    
    specific_times = []
    for _ in range(10):
        t, _ = folder_specific.monte_carlo_barrier_crossing()
        specific_times.append(t)
    
    default_mean = np.mean(default_times)
    specific_mean = np.mean(specific_times)
    
    print(f"  Default: {default_mean:.1f} ± {np.std(default_times):.1f} μs (ratio: {default_mean/exp_time_us:.2f})")
    print(f"  Specific: {specific_mean:.1f} ± {np.std(specific_times):.1f} μs (ratio: {specific_mean/exp_time_us:.2f})")
    
    # Check if we improved
    default_error = abs(np.log10(default_mean/exp_time_us))
    specific_error = abs(np.log10(specific_mean/exp_time_us))
    
    if specific_error < default_error:
        print(f"\n✅ IMPROVED by {default_error/specific_error:.1f}x (closer to experiment)")
    else:
        print(f"\n❌ Worse by {specific_error/default_error:.1f}x")


if __name__ == "__main__":
    np.random.seed(42)
    
    # Test all proteins
    proteins = [
        ('Trp-cage', 'NLYIQWLKDGGPSSGRPPPS', 4.1, 296),
        ('Villin', 'LSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF', 0.7, 300),
        ('BBA5', 'EQYTAKYKGRTFRNEKELRDFIE', 13.0, 298),
        ('WW domain', 'GSKLPPGWEKRMSRSSGRVYYFNHITNASQWERPS', 13.0, 298)
    ]
    
    for name, seq, exp_time, exp_temp in proteins:
        test_protein_final(name, seq, exp_time, exp_temp)
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print("Protein-specific parameters derived from:")
    print("1. Barrier height from structure complexity (2-3 coins)")
    print("2. P_geom from helix/sheet content")
    print("3. P_ledger from number of domains")
    print("4. All from first principles - NO fitting!") 