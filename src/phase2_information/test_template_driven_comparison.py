#!/usr/bin/env python3
"""
Compare default parameters (v2) with template-driven parameters (v3).
Tests the full protein suite to see if template analysis improves predictions.
"""

import numpy as np
from accelerated_folder_v2 import AcceleratedFolderV2
from accelerated_folder_v3 import AcceleratedFolderV3
import matplotlib.pyplot as plt

# Test proteins with experimental data
TEST_PROTEINS = [
    {
        'name': 'Trp-cage (1L2Y)',
        'sequence': 'NLYIQWLKDGGPSSGRPPPS',
        'exp_time_us': 4.1,
        'temp_K': 296,
        'reference': 'Qiu et al., 2002'
    },
    {
        'name': 'BBA5',
        'sequence': 'EQYTAKYKGRTFRNEKELRDFIE',
        'exp_time_us': 13.0,
        'temp_K': 298,
        'reference': 'Dimitriadis et al., 2004'
    },
    {
        'name': 'WW domain (FiP35)',
        'sequence': 'GSKLPPGWEKRMSRSSGRVYYFNHITNASQWERPS',
        'exp_time_us': 13.0,
        'temp_K': 298,
        'reference': 'Liu et al., 2008'
    },
    {
        'name': 'Villin HP35',
        'sequence': 'LSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF',
        'exp_time_us': 0.7,
        'temp_K': 300,
        'reference': 'Kubelka et al., 2003'
    },
    {
        'name': 'Protein G B1',
        'sequence': 'MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE',
        'exp_time_us': 170.0,
        'temp_K': 298,
        'reference': 'Park et al., 1999'
    }
]

def test_protein_comparison(protein_info, n_runs=5):
    """Test both v2 and v3 on a single protein."""
    print(f"\n{'='*70}")
    print(f"Testing {protein_info['name']}")
    print(f"{'='*70}")
    print(f"Sequence ({len(protein_info['sequence'])} residues): {protein_info['sequence']}")
    print(f"Experimental: {protein_info['exp_time_us']} μs at {protein_info['temp_K']}K")
    print(f"Reference: {protein_info['reference']}")
    
    # Test with v2 (default parameters)
    print("\n--- Testing with DEFAULT parameters (v2) ---")
    folder_v2 = AcceleratedFolderV2(
        n_residues=len(protein_info['sequence']),
        sequence=protein_info['sequence'],
        temperature=protein_info['temp_K'],
        use_protein_specific=False  # Use default parameters
    )
    
    v2_times = []
    for i in range(n_runs):
        result = folder_v2.run_accelerated()
        if result.get('mc_folding_time_us') is not None:
            v2_times.append(result['mc_folding_time_us'])
            print(f"  Run {i+1}: {v2_times[-1]:.1f} μs")
    
    # Test with v3 (template-driven parameters)
    print("\n--- Testing with TEMPLATE-DRIVEN parameters (v3) ---")
    folder_v3 = AcceleratedFolderV3(
        n_residues=len(protein_info['sequence']),
        sequence=protein_info['sequence'],
        temperature=protein_info['temp_K']
    )
    
    v3_times = []
    v3_params = None
    for i in range(n_runs):
        result = folder_v3.run_accelerated()
        if result.get('mc_folding_time_us') is not None:
            v3_times.append(result['mc_folding_time_us'])
            print(f"  Run {i+1}: {v3_times[-1]:.1f} μs")
            if v3_params is None and 'parameters' in result:
                v3_params = result['parameters']
    
    # Print parameter comparison
    if v3_params:
        print("\nTemplate-derived parameters:")
        print(f"  Barrier: {v3_params.get('barrier_eV', 0.18):.2f} eV (should always be 0.18)")
        print(f"  P_ledger: {v3_params.get('p_ledger', 1.0):.3f}")
        print(f"  P_geom: {v3_params.get('p_geom', 1.0):.3f}")
        print(f"  Path entropy: {v3_params.get('path_entropy', 1.0):.3f}")
        print(f"  Mobility anisotropy: {v3_params.get('mobility_anisotropy', 0.0):.3f}")
    
    # Calculate statistics
    results = {
        'name': protein_info['name'],
        'exp_time': protein_info['exp_time_us'],
        'v2_mean': np.mean(v2_times) if v2_times else None,
        'v2_std': np.std(v2_times) if v2_times else None,
        'v3_mean': np.mean(v3_times) if v3_times else None,
        'v3_std': np.std(v3_times) if v3_times else None,
    }
    
    if results['v2_mean'] and results['v3_mean']:
        results['v2_ratio'] = results['v2_mean'] / protein_info['exp_time_us']
        results['v3_ratio'] = results['v3_mean'] / protein_info['exp_time_us']
        results['improvement'] = results['v2_ratio'] / results['v3_ratio']
        
        print(f"\nResults:")
        print(f"  Experimental: {protein_info['exp_time_us']} μs")
        print(f"  V2 (default): {results['v2_mean']:.1f} ± {results['v2_std']:.1f} μs (ratio: {results['v2_ratio']:.2f})")
        print(f"  V3 (template): {results['v3_mean']:.1f} ± {results['v3_std']:.1f} μs (ratio: {results['v3_ratio']:.2f})")
        print(f"  Improvement: {results['improvement']:.2f}x closer to experiment")
    
    return results

def main():
    """Run comparison tests on all proteins."""
    print("="*70)
    print("TEMPLATE-DRIVEN PARAMETER COMPARISON TEST")
    print("="*70)
    print("Comparing default parameters (v2) with template-driven parameters (v3)")
    print("Goal: See if template analysis improves folding time predictions")
    
    all_results = []
    for protein in TEST_PROTEINS:
        results = test_protein_comparison(protein, n_runs=5)
        all_results.append(results)
    
    # Summary table
    print("\n" + "="*70)
    print("SUMMARY TABLE")
    print("="*70)
    print(f"{'Protein':<20} {'Exp (μs)':<10} {'V2 Ratio':<10} {'V3 Ratio':<10} {'Improved?':<10}")
    print("-"*70)
    
    improvements = []
    for r in all_results:
        if r['v2_mean'] and r['v3_mean']:
            improved = "YES" if r['improvement'] > 1.1 else "NO"
            if r['improvement'] > 1.1:
                improvements.append(r['improvement'])
            print(f"{r['name']:<20} {r['exp_time']:<10.1f} {r['v2_ratio']:<10.2f} {r['v3_ratio']:<10.2f} {improved:<10}")
    
    if improvements:
        print(f"\nAverage improvement for proteins that improved: {np.mean(improvements):.2f}x")
    
    # Create comparison plot
    plt.figure(figsize=(10, 6))
    proteins = [r['name'] for r in all_results if r['v2_mean'] and r['v3_mean']]
    v2_ratios = [r['v2_ratio'] for r in all_results if r['v2_mean'] and r['v3_mean']]
    v3_ratios = [r['v3_ratio'] for r in all_results if r['v2_mean'] and r['v3_mean']]
    
    x = np.arange(len(proteins))
    width = 0.35
    
    plt.bar(x - width/2, v2_ratios, width, label='Default (v2)', alpha=0.8)
    plt.bar(x + width/2, v3_ratios, width, label='Template-driven (v3)', alpha=0.8)
    
    plt.axhline(y=1, color='k', linestyle='--', label='Perfect prediction')
    plt.axhline(y=2, color='r', linestyle=':', alpha=0.5, label='Factor of 2')
    plt.axhline(y=0.5, color='r', linestyle=':', alpha=0.5)
    
    plt.xlabel('Protein')
    plt.ylabel('Ratio (RS/Experimental)')
    plt.title('Template-Driven Parameters: Improvement Analysis')
    plt.xticks(x, proteins, rotation=45, ha='right')
    plt.yscale('log')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('template_comparison_results.png', dpi=150)
    print(f"\nComparison plot saved to template_comparison_results.png")

if __name__ == "__main__":
    main() 