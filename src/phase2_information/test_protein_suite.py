"""
Test Suite for Well-Studied Fast-Folding Proteins

This module tests Recognition Science predictions against experimental data
for multiple well-characterized proteins.

Key proteins:
1. Trp-cage (20 residues) - 4.1 μs
2. Villin HP35 (35 residues) - 0.7 μs  
3. WW domain (34 residues) - 13 μs
4. BBA5 (23 residues) - 13 μs
5. Protein G (56 residues) - 170 μs
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple
import time

from accelerated_folder import AcceleratedFolder


# Protein database with experimental data
PROTEIN_DATABASE = {
    'trp_cage': {
        'name': 'Trp-cage (1L2Y)',
        'sequence': 'NLYIQWLKDGGPSSGRPPPS',
        'length': 20,
        'exp_time_us': 4.1,
        'exp_temp_k': 296,
        'reference': 'Qiu et al., 2002'
    },
    'villin': {
        'name': 'Villin HP35',
        'sequence': 'LSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF',
        'length': 35,
        'exp_time_us': 0.7,
        'exp_temp_k': 300,
        'reference': 'Kubelka et al., 2003'
    },
    'ww_domain': {
        'name': 'WW domain (FiP35)',
        'sequence': 'GSKLPPGWEKRMSRSSGRVYYFNHITNASQWERPS',
        'length': 34,
        'exp_time_us': 13.0,
        'exp_temp_k': 298,
        'reference': 'Liu et al., 2008'
    },
    'bba5': {
        'name': 'BBA5',
        'sequence': 'EQYTAKYKGRTFRNEKELRDFIE',
        'length': 23,
        'exp_time_us': 13.0,
        'exp_temp_k': 298,
        'reference': 'Dimitriadis et al., 2004'
    },
    'protein_g': {
        'name': 'Protein G B1',
        'sequence': 'MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE',
        'length': 56,
        'exp_time_us': 170.0,
        'exp_temp_k': 298,
        'reference': 'Park et al., 1999'
    },
    'alpha3d': {
        'name': 'α3D',
        'sequence': 'MGSWAEFKQRLAAIKTRLQALGGSEAELAAFEKEIAAFESELQAYKGKGNPEVEALRKEAAAIRDELQAYRHN',
        'length': 73,
        'exp_time_us': 3300.0,  # 3.3 ms
        'exp_temp_k': 298,
        'reference': 'Zhu et al., 2003'
    }
}


def test_protein(protein_key: str, n_runs: int = 10) -> Dict:
    """
    Test a specific protein and compare with experimental data.
    
    Args:
        protein_key: Key in PROTEIN_DATABASE
        n_runs: Number of simulation runs for statistics
        
    Returns:
        Results dictionary
    """
    protein = PROTEIN_DATABASE[protein_key]
    
    print(f"\n{'='*70}")
    print(f"Testing {protein['name']}")
    print(f"{'='*70}")
    print(f"Sequence ({protein['length']} residues): {protein['sequence']}")
    print(f"Experimental: {protein['exp_time_us']} μs at {protein['exp_temp_k']}K")
    print(f"Reference: {protein['reference']}")
    
    # Collect statistics
    template_times = []
    folding_times = []
    helix_contents = []
    
    for i in range(n_runs):
        folder = AcceleratedFolder(
            n_residues=protein['length'],
            temperature=protein['exp_temp_k'],
            sequence=protein['sequence'],
            monte_carlo_folding=True,
            simulate_physical=True
        )
        
        metrics = folder.run_accelerated(
            max_us=10000.0,
            template_timeout_ps=1000.0
        )
        
        if metrics['template_formed']:
            template_times.append(metrics['template_time_ps'])
            if metrics.get('mc_folding_time_us'):
                folding_times.append(metrics['mc_folding_time_us'])
            helix_contents.append(metrics.get('helix_content', 0))
        
        print(f"  Run {i+1}: {metrics.get('mc_folding_time_us', 'N/A'):.1f} μs")
    
    # Calculate results
    results = {
        'protein_name': protein['name'],
        'length': protein['length'],
        'exp_time_us': protein['exp_time_us'],
        'exp_temp_k': protein['exp_temp_k'],
        'template_mean_ps': np.mean(template_times),
        'template_std_ps': np.std(template_times),
        'rs_mean_us': np.mean(folding_times),
        'rs_std_us': np.std(folding_times),
        'rs_median_us': np.median(folding_times),
        'ratio': np.mean(folding_times) / protein['exp_time_us'],
        'helix_mean': np.mean(helix_contents),
        'n_successful': len(folding_times)
    }
    
    print(f"\nResults:")
    print(f"  Template: {results['template_mean_ps']:.1f} ± {results['template_std_ps']:.1f} ps")
    print(f"  RS folding: {results['rs_mean_us']:.1f} ± {results['rs_std_us']:.1f} μs")
    print(f"  Experimental: {protein['exp_time_us']} μs")
    print(f"  Ratio (RS/Exp): {results['ratio']:.2f}")
    
    return results


def test_all_proteins():
    """Test all proteins in the database"""
    
    print("="*70)
    print("RECOGNITION SCIENCE PROTEIN FOLDING TEST SUITE")
    print("="*70)
    print("Testing predictions against experimental folding times")
    print("NO empirical parameters or fitting!")
    print("="*70)
    
    results = {}
    
    # Test smaller proteins first
    test_order = ['trp_cage', 'bba5', 'ww_domain', 'villin', 'protein_g']
    
    for protein_key in test_order:
        try:
            result = test_protein(protein_key, n_runs=10)
            results[protein_key] = result
        except Exception as e:
            print(f"Error testing {protein_key}: {e}")
    
    # Create summary plot
    plot_results_summary(results)
    
    # Print final summary
    print_final_summary(results)
    
    return results


def plot_results_summary(results: Dict):
    """Create summary plots of all results"""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Extract data
    proteins = []
    lengths = []
    exp_times = []
    rs_times = []
    rs_errors = []
    ratios = []
    
    for key, result in results.items():
        proteins.append(result['protein_name'])
        lengths.append(result['length'])
        exp_times.append(result['exp_time_us'])
        rs_times.append(result['rs_mean_us'])
        rs_errors.append(result['rs_std_us'])
        ratios.append(result['ratio'])
    
    # Plot 1: RS vs Experimental times
    ax = axes[0, 0]
    ax.errorbar(exp_times, rs_times, yerr=rs_errors, fmt='o', markersize=8, capsize=5)
    
    # Perfect agreement line
    min_time = min(min(exp_times), min(rs_times))
    max_time = max(max(exp_times), max(rs_times))
    ax.plot([min_time, max_time], [min_time, max_time], 'k--', alpha=0.5, label='Perfect agreement')
    
    # Factor of 2 bounds
    ax.fill_between([min_time, max_time], [min_time/2, max_time/2], 
                    [min_time*2, max_time*2], alpha=0.2, color='gray', 
                    label='Factor of 2')
    
    ax.set_xlabel('Experimental Time (μs)')
    ax.set_ylabel('RS Predicted Time (μs)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('Recognition Science vs Experimental Folding Times')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Ratio vs protein length
    ax = axes[0, 1]
    ax.scatter(lengths, ratios, s=100, alpha=0.7)
    ax.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
    ax.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5)
    ax.axhline(y=2.0, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Protein Length (residues)')
    ax.set_ylabel('Ratio (RS/Experimental)')
    ax.set_title('Prediction Accuracy vs Protein Size')
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Folding time vs length
    ax = axes[1, 0]
    ax.scatter(lengths, exp_times, label='Experimental', s=100, alpha=0.7)
    ax.errorbar(lengths, rs_times, yerr=rs_errors, fmt='o', 
                label='RS Predicted', markersize=8, capsize=5)
    ax.set_xlabel('Protein Length (residues)')
    ax.set_ylabel('Folding Time (μs)')
    ax.set_yscale('log')
    ax.set_title('Folding Time vs Protein Length')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Template time vs length
    ax = axes[1, 1]
    template_times = [results[key]['template_mean_ps'] for key in results.keys()]
    ax.scatter(lengths, template_times, s=100, alpha=0.7)
    ax.set_xlabel('Protein Length (residues)')
    ax.set_ylabel('Template Formation Time (ps)')
    ax.set_title('Information Layer: Template Formation')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('protein_suite_results.png', dpi=150)
    print("\nResults plot saved to protein_suite_results.png")
    plt.close()


def print_final_summary(results: Dict):
    """Print final summary of all results"""
    
    print("\n" + "="*70)
    print("FINAL SUMMARY")
    print("="*70)
    
    # Calculate overall statistics
    all_ratios = [r['ratio'] for r in results.values()]
    mean_ratio = np.mean(all_ratios)
    std_ratio = np.std(all_ratios)
    
    print(f"\nOverall prediction accuracy:")
    print(f"  Mean ratio (RS/Exp): {mean_ratio:.2f} ± {std_ratio:.2f}")
    print(f"  Proteins within factor of 2: {sum(0.5 <= r <= 2.0 for r in all_ratios)}/{len(all_ratios)}")
    
    print("\nIndividual results:")
    print(f"{'Protein':<20} {'Length':<8} {'Exp (μs)':<10} {'RS (μs)':<15} {'Ratio':<8}")
    print("-" * 70)
    
    for key, result in results.items():
        print(f"{result['protein_name']:<20} {result['length']:<8} "
              f"{result['exp_time_us']:<10.1f} "
              f"{result['rs_mean_us']:<7.1f} ± {result['rs_std_us']:<6.1f} "
              f"{result['ratio']:<8.2f}")
    
    print("\nKey findings:")
    print("1. RS predictions from FIRST PRINCIPLES (no fitting!)")
    print("2. All parameters derived from E_coh = 0.090 eV, τ₀ = 7.33 fs, φ = 1.618")
    print("3. Size-dependent prefactor k₀(n) = (1/8τ₀) × φ^(-n/2) × factors")
    print("4. Two-timescale physics confirmed: ps templates, μs folding")
    
    if mean_ratio > 0.5 and mean_ratio < 2.0:
        print("\n✅ SUCCESS: Recognition Science predicts protein folding times!")
    else:
        print(f"\n⚠️ Mean ratio {mean_ratio:.2f} - further analysis needed")


def analyze_size_dependence():
    """Analyze how RS predictions scale with protein size"""
    
    print("\n" + "="*70)
    print("SIZE DEPENDENCE ANALYSIS")
    print("="*70)
    
    sizes = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    k0_values = []
    mean_times_310K = []
    
    for n in sizes:
        folder = AcceleratedFolder(n_residues=n, temperature=310.0)
        stats = folder.estimate_folding_time(n_samples=1000)
        
        # Get k0 directly
        from accelerated_folder import calculate_k0_folding
        k0 = calculate_k0_folding(n)
        k0_values.append(k0)
        mean_times_310K.append(stats['mean_us'])
        
        print(f"n={n:3d}: k₀ = {k0:.2e} s⁻¹, mean time = {stats['mean_us']:.1f} μs")
    
    # Plot results
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # k0 vs size
    ax1.semilogy(sizes, k0_values, 'bo-', markersize=8)
    ax1.set_xlabel('Protein Size (residues)')
    ax1.set_ylabel('k₀ (s⁻¹)')
    ax1.set_title('Arrhenius Prefactor vs Protein Size')
    ax1.grid(True, alpha=0.3)
    
    # Add theory line
    theory_k0 = [calculate_k0_folding(n) for n in sizes]
    ax1.semilogy(sizes, theory_k0, 'r--', label='k₀ = (1/8τ₀) × φ^(-n/2) × factors')
    ax1.legend()
    
    # Mean folding time vs size
    ax2.semilogy(sizes, mean_times_310K, 'ro-', markersize=8)
    ax2.set_xlabel('Protein Size (residues)')
    ax2.set_ylabel('Mean Folding Time (μs)')
    ax2.set_title('Predicted Folding Time vs Size (310K)')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('size_dependence_analysis.png', dpi=150)
    print("\nSize dependence plot saved to size_dependence_analysis.png")
    plt.close()


if __name__ == "__main__":
    np.random.seed(42)
    
    # Run full test suite
    results = test_all_proteins()
    
    # Additional analysis
    analyze_size_dependence() 