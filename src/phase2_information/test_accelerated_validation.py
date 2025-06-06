"""
Validation Test: Compare Accelerated vs Normal Simulation

This script validates that the accelerated folder produces
statistically consistent results with the normal simulation.
"""

import numpy as np
from typing import Dict, List
import time

from enhanced_three_layer_folder import EnhancedThreeLayerFolder
from accelerated_folder import AcceleratedFolder


def run_normal_simulation(n_residues: int, sequence: str, max_us: float = 10.0) -> Dict:
    """Run normal (non-accelerated) simulation"""
    print(f"\nRunning NORMAL simulation for {n_residues} residues...")
    
    folder = EnhancedThreeLayerFolder(
        n_residues=n_residues,
        temperature=310.0,
        sequence=sequence
    )
    
    start_time = time.time()
    max_ticks = int(max_us * 1e-6 / 7.33e-15)
    
    # Track key events
    template_time = None
    folding_time = None
    
    for tick in range(max_ticks):
        folder.step()
        
        # Check for template
        if folder.template_ready and template_time is None:
            template_time = folder.template_completion_tick * 7.33e-15 * 1e12  # ps
            print(f"  Template formed at {template_time:.1f} ps")
        
        # Check for folding
        if folder.folding_initiated and folding_time is None:
            folding_time = folder.tick * 7.33e-15 * 1e6  # μs
            print(f"  Folding initiated at {folding_time:.2f} μs")
            break  # Stop here for comparison
        
        # Progress updates
        if tick % 100000 == 0 and tick > 0:
            current_us = tick * 7.33e-15 * 1e6
            print(f"  Progress: {current_us:.2f} μs")
    
    wall_time = time.time() - start_time
    
    return {
        'template_time_ps': template_time,
        'folding_time_us': folding_time,
        'wall_time_s': wall_time,
        'recognition_events': len(folder.recognition_events),
        'final_coherence': folder.phase_field.get_folding_progress()['coherence']
    }


def run_accelerated_simulation(n_residues: int, sequence: str) -> Dict:
    """Run accelerated simulation"""
    print(f"\nRunning ACCELERATED simulation for {n_residues} residues...")
    
    folder = AcceleratedFolder(
        n_residues=n_residues,
        temperature=310.0,
        sequence=sequence,
        monte_carlo_folding=True,
        simulate_physical=False  # Just test barrier crossing
    )
    
    metrics = folder.run_accelerated(
        max_us=10000.0,
        template_timeout_ps=1000.0
    )
    
    return {
        'template_time_ps': metrics['template_time_ps'],
        'folding_time_us': metrics['mc_folding_time_us'],
        'wall_time_s': metrics['wall_time_s'],
        'recognition_events': metrics['template_recognitions'],
        'final_coherence': metrics.get('phase_coherence', 0)
    }


def compare_results(normal: Dict, accelerated: Dict):
    """Compare normal vs accelerated results"""
    print("\n" + "="*50)
    print("COMPARISON: Normal vs Accelerated")
    print("="*50)
    
    print(f"\nTemplate Formation:")
    print(f"  Normal:      {normal['template_time_ps']:.1f} ps")
    print(f"  Accelerated: {accelerated['template_time_ps']:.1f} ps")
    print(f"  Difference:  {abs(normal['template_time_ps'] - accelerated['template_time_ps']):.1f} ps")
    
    print(f"\nFolding Time:")
    if normal['folding_time_us'] and accelerated['folding_time_us']:
        print(f"  Normal:      {normal['folding_time_us']:.2f} μs")
        print(f"  Accelerated: {accelerated['folding_time_us']:.2f} μs")
        print(f"  Ratio:       {accelerated['folding_time_us'] / normal['folding_time_us']:.2f}")
    else:
        print("  One or both simulations didn't reach folding")
    
    print(f"\nComputational Performance:")
    print(f"  Normal:      {normal['wall_time_s']:.1f} s")
    print(f"  Accelerated: {accelerated['wall_time_s']:.1f} s")
    print(f"  Speedup:     {normal['wall_time_s'] / accelerated['wall_time_s']:.0f}x")
    
    print(f"\nRecognition Events (template phase):")
    print(f"  Normal:      {normal['recognition_events']}")
    print(f"  Accelerated: {accelerated['recognition_events']}")


def test_statistical_consistency():
    """Test that MC folding times match expected distribution"""
    print("\n" + "="*50)
    print("STATISTICAL VALIDATION")
    print("="*50)
    
    n_residues = 10
    folder = AcceleratedFolder(n_residues=n_residues, temperature=310.0)
    
    # Get theoretical rate
    k_fold = 3.16e6 * np.exp(-0.18 / folder.kT)
    expected_mean = 1e6 / k_fold  # μs
    
    # Sample many folding times
    n_samples = 1000
    folding_times = []
    for _ in range(n_samples):
        time_us, _ = folder.monte_carlo_barrier_crossing()
        folding_times.append(time_us)
    
    # Statistics
    mean_time = np.mean(folding_times)
    std_time = np.std(folding_times)
    
    print(f"\nMonte Carlo Folding Times (n={n_samples}):")
    print(f"  Expected mean: {expected_mean:.1f} μs")
    print(f"  Actual mean:   {mean_time:.1f} μs")
    print(f"  Std dev:       {std_time:.1f} μs")
    print(f"  Error:         {abs(mean_time - expected_mean) / expected_mean * 100:.1f}%")
    
    # Check if distribution is exponential
    # For exponential distribution, mean ≈ std dev
    print(f"\nDistribution check:")
    print(f"  Mean/StdDev ratio: {mean_time/std_time:.2f} (expect ~1.0 for exponential)")


def main():
    """Run validation tests"""
    print("="*60)
    print("ACCELERATED FOLDER VALIDATION")
    print("="*60)
    
    # Test 1: Small peptide comparison
    print("\nTest 1: Small peptide (5 residues)")
    sequence = "AAAAA"
    
    normal_result = run_normal_simulation(5, sequence, max_us=1.0)
    accel_result = run_accelerated_simulation(5, sequence)
    compare_results(normal_result, accel_result)
    
    # Test 2: Statistical validation
    test_statistical_consistency()
    
    # Test 3: Check template formation consistency
    print("\n" + "="*50)
    print("TEMPLATE FORMATION CONSISTENCY")
    print("="*50)
    
    # Run same sequence multiple times
    n_runs = 5
    template_times = []
    
    for i in range(n_runs):
        folder = AcceleratedFolder(n_residues=5, temperature=310.0, sequence="AAAAA")
        metrics = folder.run_until_template(max_ps=1000.0)
        if metrics['template_formed']:
            template_times.append(metrics['template_time_ps'])
    
    if template_times:
        print(f"\nTemplate times over {n_runs} runs:")
        for i, t in enumerate(template_times):
            print(f"  Run {i+1}: {t:.1f} ps")
        print(f"  Mean: {np.mean(template_times):.1f} ps")
        print(f"  Std:  {np.std(template_times):.1f} ps")
    
    print("\n" + "="*60)
    print("VALIDATION COMPLETE")
    print("="*60)
    print("\nConclusion: The accelerated folder should produce")
    print("statistically consistent results while being much faster!")


if __name__ == "__main__":
    np.random.seed(42)
    main() 