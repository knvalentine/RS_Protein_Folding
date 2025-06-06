"""
Test Larger Proteins with Accelerated Folder

This script demonstrates how we can now test 15-30 residue proteins
efficiently using Monte Carlo barrier crossing.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict
import time

from accelerated_folder import AcceleratedFolder


def test_protein_sequences():
    """Test a range of protein sequences from small to medium"""
    
    sequences = [
        # Small peptides (control)
        ("AAAAA", "Penta-alanine", "helix"),
        ("VVVVV", "Penta-valine", "sheet"),
        
        # 10-residue peptides
        ("AAAAAEAAAA", "Helix-10", "helix"),
        ("VVVVVVVVVV", "Sheet-10", "sheet"),
        ("GGGGGGGGGG", "Flexible-10", "coil"),
        
        # 15-residue peptides
        ("AAAAAAAAAAEAAAA", "Helix-15", "helix"),
        ("VVVVVVVVVVVVVVV", "Sheet-15", "sheet"),
        ("AEAEAEAEAEAEAEA", "Mixed-15", "mixed"),
        
        # 20-residue peptides (real test)
        ("AAAAAAAAAAAAAAAAAAAA", "Helix-20", "helix"),
        ("VVVVVVVVVVVVVVVVVVVV", "Sheet-20", "sheet"),
        ("AEAEAEAEAEAEAEAEAEAE", "Mixed-20", "mixed"),
        
        # 25-residue peptides
        ("A" * 25, "Helix-25", "helix"),
        ("V" * 25, "Sheet-25", "sheet"),
        
        # 30-residue peptide
        ("A" * 30, "Helix-30", "helix"),
    ]
    
    print("="*70)
    print("TESTING LARGER PROTEINS WITH ACCELERATED FOLDER")
    print("="*70)
    
    results = []
    
    for sequence, name, expected_structure in sequences:
        print(f"\n{'='*50}")
        print(f"Testing: {name}")
        print(f"Sequence: {sequence} (n={len(sequence)})")
        print(f"Expected: {expected_structure}")
        print(f"{'='*50}")
        
        # Create accelerated folder
        folder = AcceleratedFolder(
            n_residues=len(sequence),
            temperature=310.0,
            sequence=sequence,
            monte_carlo_folding=True,
            simulate_physical=True
        )
        
        # Run accelerated simulation
        metrics = folder.run_accelerated(
            max_us=10000.0,  # 10 ms max
            template_timeout_ps=1000.0  # 1 ns for template
        )
        
        # Store results
        results.append({
            'name': name,
            'sequence': sequence,
            'n_residues': len(sequence),
            'expected': expected_structure,
            **metrics
        })
        
        # Print summary
        print(f"\nResults:")
        print(f"  Template time: {metrics.get('template_time_ps', 'N/A'):.1f} ps")
        print(f"  MC folding time: {metrics.get('mc_folding_time_us', 'N/A'):.1f} μs")
        print(f"  Helix content: {metrics.get('helix_content', 0)*100:.1f}%")
        print(f"  Sheet content: {metrics.get('sheet_content', 0)*100:.1f}%")
        print(f"  Wall time: {metrics.get('wall_time_s', 0):.1f} s")
        print(f"  Acceleration: {metrics.get('acceleration_factor', 0):.0f}x")
    
    return results


def analyze_scaling(results: List[Dict]):
    """Analyze how metrics scale with protein size"""
    
    print("\n" + "="*70)
    print("SCALING ANALYSIS")
    print("="*70)
    
    # Extract data
    n_residues = [r['n_residues'] for r in results if r['template_formed']]
    template_times = [r['template_time_ps'] for r in results if r['template_formed']]
    folding_times = [r['mc_folding_time_us'] for r in results if r['mc_folding_time_us']]
    wall_times = [r['wall_time_s'] for r in results if r['template_formed']]
    
    # Print scaling summary
    print(f"\n{'N':<5} {'Template(ps)':<12} {'Folding(μs)':<12} {'Wall(s)':<10}")
    print("-"*40)
    
    for r in results:
        if r['template_formed']:
            print(f"{r['n_residues']:<5} "
                  f"{r['template_time_ps']:<12.1f} "
                  f"{r.get('mc_folding_time_us', 0):<12.1f} "
                  f"{r['wall_time_s']:<10.1f}")
    
    # Plot scaling
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # Template time vs size
    ax1.scatter(n_residues, template_times, color='blue', s=100, alpha=0.7)
    ax1.set_xlabel('Number of Residues')
    ax1.set_ylabel('Template Formation Time (ps)')
    ax1.set_title('Information Layer Scaling')
    ax1.grid(True, alpha=0.3)
    
    # Folding time vs size
    if folding_times:
        ax2.scatter([r['n_residues'] for r in results if r.get('mc_folding_time_us')],
                   folding_times, color='red', s=100, alpha=0.7)
        ax2.set_xlabel('Number of Residues')
        ax2.set_ylabel('Folding Time (μs)')
        ax2.set_title('Physical Layer Scaling')
        ax2.set_yscale('log')
        ax2.grid(True, alpha=0.3)
    
    # Wall time vs size
    ax3.scatter(n_residues, wall_times, color='green', s=100, alpha=0.7)
    ax3.set_xlabel('Number of Residues')
    ax3.set_ylabel('Computation Time (s)')
    ax3.set_title('Computational Efficiency')
    ax3.grid(True, alpha=0.3)
    
    # Acceleration factor
    accelerations = [r.get('acceleration_factor', 0) for r in results if r['template_formed']]
    ax4.scatter(n_residues, accelerations, color='purple', s=100, alpha=0.7)
    ax4.set_xlabel('Number of Residues')
    ax4.set_ylabel('Acceleration Factor')
    ax4.set_title('Speedup vs Direct Simulation')
    ax4.set_yscale('log')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('accelerated_scaling_analysis.png', dpi=150)
    print("\nScaling plots saved to accelerated_scaling_analysis.png")
    plt.close()


def test_folding_statistics():
    """Test Monte Carlo folding time statistics"""
    
    print("\n" + "="*70)
    print("MONTE CARLO FOLDING TIME STATISTICS")
    print("="*70)
    
    test_sizes = [10, 20, 30, 40]
    
    for n in test_sizes:
        folder = AcceleratedFolder(n_residues=n, temperature=310.0)
        stats = folder.estimate_folding_time(n_samples=10000)
        
        print(f"\nN = {n} residues:")
        print(f"  Rate constant: {stats['rate_constant']:.2e} s^-1")
        print(f"  Mean folding time: {stats['mean_us']:.1f} μs")
        print(f"  Median: {stats['median_us']:.1f} μs")
        print(f"  Std dev: {stats['std_us']:.1f} μs")
        print(f"  Range: {stats['min_us']:.1f} - {stats['max_us']:.1f} μs")


def compare_with_theory():
    """Compare results with RS theoretical predictions"""
    
    print("\n" + "="*70)
    print("COMPARISON WITH RS THEORY")
    print("="*70)
    
    # Theoretical predictions
    print("\nTheoretical RS Predictions:")
    print("1. Template formation: ~65 ps (independent of size)")
    print("2. Folding barrier: 0.18 eV = 2 × E_coh")
    print("3. Prefactor scales as φ^(-n/2)")
    print("4. Two distinct timescales: ps vs μs")
    
    # Test specific size
    n = 40  # Where prefactor ≈ 3.16×10^6
    folder = AcceleratedFolder(n_residues=n, temperature=310.0)
    stats = folder.estimate_folding_time()
    
    print(f"\nFor n={n} residues:")
    print(f"  Expected k₀: 3.16×10^6 s^-1")
    print(f"  Actual k: {stats['rate_constant']:.2e} s^-1")
    print(f"  Expected mean time: ~280 μs")
    print(f"  Actual mean time: {stats['mean_us']:.1f} μs")


if __name__ == "__main__":
    # Set random seed
    np.random.seed(42)
    
    # Run tests
    results = test_protein_sequences()
    
    # Analysis
    analyze_scaling(results)
    test_folding_statistics()
    compare_with_theory()
    
    print("\n" + "="*70)
    print("ACCELERATED TESTING COMPLETE")
    print("="*70)
    print("\nKey Achievement: We can now test 20-30 residue proteins efficiently!")
    print("Monte Carlo barrier crossing enables realistic protein size studies.") 