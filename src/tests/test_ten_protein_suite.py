#!/usr/bin/env python3
"""
Test the full 10-protein suite with template-driven parameters.
This establishes a baseline for the current implementation.
"""

import numpy as np
from accelerated_folder_v3 import AcceleratedFolderV3
import matplotlib.pyplot as plt
import time

# Full 10-protein test suite
TEST_PROTEINS = [
    {
        'name': 'Trp-cage',
        'pdb': '1L2Y',
        'sequence': 'NLYIQWLKDGGPSSGRPPPS',
        'exp_time_us': 4.1,
        'temp_K': 296,
        'reference': 'Qiu et al., 2002'
    },
    {
        'name': 'WW domain',
        'pdb': 'FiP35',
        'sequence': 'GSKLPPGWEKRMSRSSGRVYYFNHITNASQWERPS',
        'exp_time_us': 13.0,
        'temp_K': 298,
        'reference': 'Liu et al., 2008'
    },
    {
        'name': 'Villin HP35',
        'pdb': '1VII',
        'sequence': 'LSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF',
        'exp_time_us': 0.7,
        'temp_K': 300,
        'reference': 'Kubelka et al., 2003'
    },
    {
        'name': 'BBA5',
        'pdb': '1FME',
        'sequence': 'EQYTAKYKGRTFRNEKELRDFIE',
        'exp_time_us': 13.0,
        'temp_K': 298,
        'reference': 'Dimitriadis et al., 2004'
    },
    {
        'name': 'Protein G B1',
        'pdb': '1GB1',
        'sequence': 'MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE',
        'exp_time_us': 170.0,
        'temp_K': 298,
        'reference': 'Park et al., 1999'
    },
    {
        'name': 'Pin1 WW',
        'pdb': '1PIN',
        'sequence': 'KLPPGWEKRMSRSSGRVYYFNHITNASQFERPSG',
        'exp_time_us': 18.0,
        'temp_K': 298,
        'reference': 'Jäger et al., 2001'
    },
    {
        'name': 'Protein L',
        'pdb': '1HZ6',
        'sequence': 'MEEVTIKANLIFANGSTQTAEFKGTFEKATSEAYAYADTLKKDNGEWTVDVADKGYTLNIKFAG',
        'exp_time_us': 18.0,
        'temp_K': 298,
        'reference': 'Scalley et al., 1997'
    },
    {
        'name': 'Ubiquitin',
        'pdb': '1UBQ',
        'sequence': 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',
        'exp_time_us': 1000.0,  # Estimate for 76-residue
        'temp_K': 298,
        'reference': 'Went & Jackson, 2005'
    },
    {
        'name': 'Cold-shock',
        'pdb': '1CSP',
        'sequence': 'MLEGKVKWFNSEKGFGFIEVEGQDDVFVHFSAIQGEGFKTLEEGQAVSFEIVEGNRGPQAANVTKEA',
        'exp_time_us': 50.0,
        'temp_K': 298,
        'reference': 'Schindler et al., 1995'
    },
    {
        'name': 'λ-repressor',
        'pdb': '1LMB',
        'sequence': 'MQEQLEDARRLKAIYEKKKNELGLSQESVADKMGMGQSGVGALFNGINALNAYNAALLAKILKVSVEEFSPSIAREIYEMYEAVSMQPSLRSEYEYPVFSHVQAGMFSPELRTFTKGDAERWVSTTKKASDSAFWLEVEGNSMTAPTGSKPSFPDGMLILVDPEQAVEPGDFCIARLGGDEFTFKKLIRDSGQVFLQPLNPQYPMIPCNESCSVVGKVIASQWPEETFG',
        'exp_time_us': 5000.0,  # Estimate for 80-residue fragment
        'temp_K': 298,
        'reference': 'Burton et al., 1997'
    }
]

def test_protein(protein_info, n_runs=3):
    """Test a single protein with template-driven parameters."""
    print(f"\n{'='*70}")
    print(f"{protein_info['name']} ({protein_info['pdb']})")
    print(f"{'='*70}")
    print(f"Sequence ({len(protein_info['sequence'])} residues): {protein_info['sequence'][:50]}...")
    print(f"Experimental: {protein_info['exp_time_us']} μs at {protein_info['temp_K']}K")
    print(f"Reference: {protein_info['reference']}")
    
    folder = AcceleratedFolderV3(
        n_residues=len(protein_info['sequence']),
        sequence=protein_info['sequence'],
        temperature=protein_info['temp_K'],
        use_template_params=True
    )
    
    times = []
    params_list = []
    
    for i in range(n_runs):
        start = time.time()
        result = folder.run_accelerated(max_us=100000.0)  # 100 ms max
        elapsed = time.time() - start
        
        if result.get('mc_folding_time_us') is not None:
            total_time = result['total_time_us']
            times.append(total_time)
            
            if i == 0:  # Print parameters once
                params = result.get('parameters_used', {})
                print(f"\nTemplate-derived parameters:")
                print(f"  Voxels: {params.get('n_voxels', 'N/A')}")
                print(f"  Components: {params.get('n_components', 'N/A')}")
                print(f"  P_ledger: {params.get('p_ledger', 'N/A'):.3f}")
                print(f"  P_geom: {params.get('p_geom', 'N/A'):.4f}")
                print(f"  Path entropy: {params.get('path_entropy', 'N/A'):.4f}")
                print(f"  Mobility anisotropy: {params.get('mobility_anisotropy', 'N/A'):.4f}")
            
            print(f"  Run {i+1}: {total_time:.1f} μs (wall time: {elapsed:.1f}s)")
    
    if times:
        mean_time = np.mean(times)
        std_time = np.std(times)
        ratio = mean_time / protein_info['exp_time_us']
        
        print(f"\nResults:")
        print(f"  RS prediction: {mean_time:.1f} ± {std_time:.1f} μs")
        print(f"  Experimental: {protein_info['exp_time_us']} μs")
        print(f"  Ratio (RS/Exp): {ratio:.2f}")
        
        if 0.5 <= ratio <= 2.0:
            print("  ✅ Within factor of 2!")
        elif 0.1 <= ratio <= 10.0:
            print("  ⚠️  Within order of magnitude")
        else:
            print(f"  ❌ Off by {max(ratio, 1/ratio):.0f}x")
        
        return {
            'name': protein_info['name'],
            'length': len(protein_info['sequence']),
            'exp_time': protein_info['exp_time_us'],
            'rs_mean': mean_time,
            'rs_std': std_time,
            'ratio': ratio,
            'success': True
        }
    else:
        print("  ❌ No successful runs")
        return {
            'name': protein_info['name'],
            'length': len(protein_info['sequence']),
            'exp_time': protein_info['exp_time_us'],
            'success': False
        }

def main():
    """Run tests on all 10 proteins."""
    print("="*70)
    print("10-PROTEIN RECOGNITION SCIENCE TEST SUITE")
    print("="*70)
    print("Testing template-driven parameter extraction")
    print("All parameters from first principles - NO fitting!")
    
    np.random.seed(42)
    
    results = []
    for protein in TEST_PROTEINS:
        result = test_protein(protein, n_runs=3)
        results.append(result)
    
    # Summary statistics
    successful = [r for r in results if r.get('success', False)]
    
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    # Table of results
    print(f"\n{'Protein':<15} {'Length':<8} {'Exp (μs)':<12} {'RS (μs)':<12} {'Ratio':<8} {'Status':<10}")
    print("-"*70)
    
    within_2x = 0
    within_10x = 0
    
    for r in successful:
        ratio = r['ratio']
        if 0.5 <= ratio <= 2.0:
            status = "✅ <2x"
            within_2x += 1
        elif 0.1 <= ratio <= 10.0:
            status = "⚠️  <10x"
            within_10x += 1
        else:
            status = "❌ >10x"
        
        print(f"{r['name']:<15} {r['length']:<8} {r['exp_time']:<12.1f} "
              f"{r['rs_mean']:<12.1f} {ratio:<8.2f} {status:<10}")
    
    print(f"\nSuccess rate: {len(successful)}/{len(results)} proteins")
    print(f"Within factor of 2: {within_2x}/{len(successful)}")
    print(f"Within order of magnitude: {within_10x + within_2x}/{len(successful)}")
    
    # Create visualization
    plt.figure(figsize=(12, 8))
    
    # Scatter plot
    lengths = [r['length'] for r in successful]
    exp_times = [r['exp_time'] for r in successful]
    rs_times = [r['rs_mean'] for r in successful]
    
    plt.scatter(exp_times, rs_times, s=100, alpha=0.7)
    
    # Perfect prediction line
    min_time = min(min(exp_times), min(rs_times))
    max_time = max(max(exp_times), max(rs_times))
    plt.plot([min_time, max_time], [min_time, max_time], 'k--', label='Perfect prediction')
    plt.fill_between([min_time, max_time], [min_time/2, max_time/2], 
                     [min_time*2, max_time*2], alpha=0.2, color='green', 
                     label='Factor of 2')
    
    # Labels
    for r in successful:
        plt.annotate(r['name'], (r['exp_time'], r['rs_mean']), 
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    plt.xlabel('Experimental folding time (μs)')
    plt.ylabel('RS prediction (μs)')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.title('Recognition Science Protein Folding Predictions\n(Template-driven parameters, NO fitting)')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('ten_protein_results.png', dpi=150)
    print(f"\nResults plot saved to ten_protein_results.png")

if __name__ == "__main__":
    main() 