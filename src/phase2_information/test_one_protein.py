"""
Test one protein at a time for quick results
"""

import numpy as np
import sys
from accelerated_folder import AcceleratedFolder

# Quick test with just 3 runs per protein
def test_single_protein(protein_name: str = 'villin'):
    """Test a single protein with minimal runs"""
    
    proteins = {
        'villin': {
            'name': 'Villin HP35',
            'sequence': 'LSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF',
            'length': 35,
            'exp_time_us': 0.7,
            'exp_temp_k': 300,
            'reference': 'Kubelka et al., 2003'
        },
        'bba5': {
            'name': 'BBA5',
            'sequence': 'EQYTAKYKGRTFRNEKELRDFIE',
            'length': 23,
            'exp_time_us': 13.0,
            'exp_temp_k': 298,
            'reference': 'Dimitriadis et al., 2004'
        },
        'ww_domain': {
            'name': 'WW domain (FiP35)',
            'sequence': 'GSKLPPGWEKRMSRSSGRVYYFNHITNASQWERPS',
            'length': 34,
            'exp_time_us': 13.0,
            'exp_temp_k': 298,
            'reference': 'Liu et al., 2008'
        }
    }
    
    if protein_name not in proteins:
        print(f"Unknown protein: {protein_name}")
        print(f"Available: {list(proteins.keys())}")
        return
    
    protein = proteins[protein_name]
    
    print(f"\n{'='*70}")
    print(f"Testing {protein['name']}")
    print(f"{'='*70}")
    print(f"Length: {protein['length']} residues")
    print(f"Experimental: {protein['exp_time_us']} μs at {protein['exp_temp_k']}K")
    print(f"Reference: {protein['reference']}")
    print("-" * 50)
    
    folding_times = []
    template_times = []
    
    # Just 3 runs for quick test
    for i in range(3):
        print(f"\nRun {i+1}/3:")
        
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
                print(f"  Folding time: {metrics['mc_folding_time_us']:.1f} μs")
    
    # Results
    print(f"\n{'='*50}")
    print("RESULTS:")
    print(f"{'='*50}")
    print(f"Template formation: {np.mean(template_times):.1f} ± {np.std(template_times):.1f} ps")
    print(f"RS folding time: {np.mean(folding_times):.1f} ± {np.std(folding_times):.1f} μs")
    print(f"Experimental: {protein['exp_time_us']} μs")
    print(f"Ratio (RS/Exp): {np.mean(folding_times)/protein['exp_time_us']:.2f}")
    
    if 0.5 <= np.mean(folding_times)/protein['exp_time_us'] <= 2.0:
        print("\n✅ Within factor of 2!")
    else:
        print(f"\n⚠️ Off by factor of {np.mean(folding_times)/protein['exp_time_us']:.1f}")


if __name__ == "__main__":
    # Check command line argument
    if len(sys.argv) > 1:
        protein = sys.argv[1]
    else:
        protein = 'villin'  # Default to Villin
    
    np.random.seed(42)
    test_single_protein(protein) 