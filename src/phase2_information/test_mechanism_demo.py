"""
Mechanism Demonstration Test

This script demonstrates the enhanced folding mechanisms without waiting
for realistic microsecond timescales. It:

1. Shows information template formation
2. Demonstrates torsion angle evolution
3. Shows secondary structure detection
4. Uses artificially lowered barriers for quick testing
5. Analyzes IR photon emissions

This is for DEMONSTRATION, not realistic folding simulation.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from phase2_information.enhanced_three_layer_folder import EnhancedThreeLayerFolder


def demonstrate_mechanisms(sequence: str = "AAAAAAAA", 
                         barrier_factor: float = 0.1) -> Dict:
    """
    Demonstrate folding mechanisms with artificially lowered barrier.
    
    Args:
        sequence: Test sequence (default: 8 alanines)
        barrier_factor: Factor to reduce barrier by (0.1 = 10x easier)
        
    Returns:
        Dictionary of results and metrics
    """
    print(f"\n{'='*60}")
    print("MECHANISM DEMONSTRATION")
    print(f"Sequence: {sequence} ({len(sequence)} residues)")
    print(f"Barrier reduction: {barrier_factor}x (demo mode)")
    print(f"{'='*60}")
    
    # Create folder
    folder = EnhancedThreeLayerFolder(
        n_residues=len(sequence),
        temperature=310.0,
        sequence=sequence
    )
    
    # Tracking
    results = {
        'sequence': sequence,
        'n_residues': len(sequence),
        'time_points': [],
        'phase_coherence': [],
        'recognition_events': [],
        'helix_content': [],
        'sheet_content': [],
        'native_contacts': [],
        'torsion_costs': [],
        'ir_photons': []
    }
    
    # Phase 1: Template formation (detailed tracking)
    print("\nPhase 1: Information Template Formation")
    print("-" * 40)
    
    template_ticks = 0
    while not folder.template_ready and template_ticks < 50000:
        metrics = folder.step()
        
        # Record every 100 ticks
        if folder.tick % 100 == 0:
            detailed = folder.get_detailed_metrics()
            results['time_points'].append(detailed['time_ps'])
            results['phase_coherence'].append(detailed['phase_coherence'])
            results['recognition_events'].append(len(folder.recognition_events))
            results['helix_content'].append(detailed['helix_content'])
            results['sheet_content'].append(detailed['sheet_content'])
            results['native_contacts'].append(detailed['native_contacts'])
            results['torsion_costs'].append(detailed['avg_torsion_cost'])
            results['ir_photons'].append(len(folder.ir_analyzer.photon_events))
            
        # Progress reports
        if folder.tick % 5000 == 0 and folder.tick > 0:
            print(f"  {folder.tick * 7.33e-15 * 1e12:.1f} ps: "
                  f"{len(folder.recognition_events)} events, "
                  f"coherence={metrics['phase_coherence']:.3f}")
        
        template_ticks += 1
    
    if folder.template_ready:
        print(f"\n✓ Template complete at {folder.template_completion_tick * 7.33e-15 * 1e12:.1f} ps")
        print(f"  Total recognition events: {len(folder.recognition_events)}")
        print(f"  Phase coherence: {metrics['phase_coherence']:.3f}")
    
    # Phase 2: Artificial barrier crossing (for demonstration)
    print("\nPhase 2: Barrier Crossing (artificially accelerated)")
    print("-" * 40)
    
    # Override the barrier height for demonstration
    original_barrier = 0.18
    demo_barrier = original_barrier * barrier_factor
    
    # Force barrier crossing after some attempts
    max_attempts = 100
    for attempt in range(max_attempts):
        folder.tick += 1000  # Jump forward
        
        # Check with reduced barrier
        k_fold = 3.16e6 * np.exp(-demo_barrier / folder.kT)
        dt = 1000 * 7.33e-15
        p_initiate = 1 - np.exp(-k_fold * dt)
        
        if np.random.random() < p_initiate or attempt > 50:  # Force after 50 attempts
            folder.folding_initiated = True
            folder.barrier_crossing_attempts = attempt + 1
            print(f"✓ Barrier crossed after {attempt + 1} attempts")
            print(f"  Time: {folder.tick * 7.33e-15 * 1e6:.3f} μs")
            break
    
    # Phase 3: Physical folding (accelerated)
    print("\nPhase 3: Physical Folding Dynamics")
    print("-" * 40)
    
    if folder.folding_initiated:
        # Run folding for a reasonable time
        folding_steps = 10000
        
        for step in range(folding_steps):
            folder.step()
            
            # Force some torsion updates to see structure formation
            if step % 100 == 0:
                folder._update_torsion_angles()
                folder._detect_secondary_structures()
            
            # Record metrics
            if step % 500 == 0:
                detailed = folder.get_detailed_metrics()
                results['time_points'].append(detailed['time_ps'])
                results['phase_coherence'].append(detailed['phase_coherence'])
                results['recognition_events'].append(len(folder.recognition_events))
                results['helix_content'].append(detailed['helix_content'])
                results['sheet_content'].append(detailed['sheet_content'])
                results['native_contacts'].append(detailed['native_contacts'])
                results['torsion_costs'].append(detailed['avg_torsion_cost'])
                results['ir_photons'].append(len(folder.ir_analyzer.photon_events))
                
                if step % 2000 == 0 and step > 0:
                    print(f"  Step {step}: Helix={detailed['helix_content']:.1%}, "
                          f"Contacts={detailed['native_contacts']}, "
                          f"Progress={detailed['folding_progress']:.1%}")
    
    # Final analysis
    print("\nFinal Structure Analysis")
    print("-" * 40)
    
    final_metrics = folder.get_detailed_metrics()
    print(f"Helix content: {final_metrics['helix_content']:.1%}")
    print(f"Sheet content: {final_metrics['sheet_content']:.1%}")
    print(f"Native contacts: {final_metrics['native_contacts']}")
    print(f"Average torsion cost: {final_metrics['avg_torsion_cost']:.2f} E_coh")
    print(f"Voxel transitions: {final_metrics['voxel_transitions']}")
    
    # Secondary structures
    if folder.ss_templates:
        print("\nSecondary structures formed:")
        for ss_id, residues in folder.ss_templates.items():
            print(f"  {ss_id}: residues {min(residues)}-{max(residues)}")
    
    # IR photon analysis
    print("\nIR Photon Emission Summary")
    print("-" * 40)
    folder.ir_analyzer.plot_emission_analysis(f"{sequence}_demo")
    print(folder.ir_analyzer.generate_report(f"{sequence} Demo"))
    
    # Store folder for further analysis
    results['folder'] = folder
    
    return results


def plot_mechanism_evolution(results: Dict):
    """Plot the evolution of various mechanisms over time"""
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Recognition Science Folding Mechanisms', fontsize=16)
    
    time_ps = results['time_points']
    time_ns = [t/1000 for t in time_ps]  # Convert to ns
    
    # Phase coherence
    axes[0, 0].plot(time_ns, results['phase_coherence'], 'b-', linewidth=2)
    axes[0, 0].set_xlabel('Time (ns)')
    axes[0, 0].set_ylabel('Phase Coherence')
    axes[0, 0].set_title('Information Organization')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].axhline(y=0.8, color='r', linestyle='--', alpha=0.5, label='Template threshold')
    axes[0, 0].legend()
    
    # Recognition events
    axes[0, 1].plot(time_ns, results['recognition_events'], 'g-', linewidth=2)
    axes[0, 1].set_xlabel('Time (ns)')
    axes[0, 1].set_ylabel('Cumulative Events')
    axes[0, 1].set_title('Recognition Events')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Secondary structure
    axes[0, 2].plot(time_ns, [h*100 for h in results['helix_content']], 'b-', 
                    linewidth=2, label='Helix')
    axes[0, 2].plot(time_ns, [s*100 for s in results['sheet_content']], 'r-', 
                    linewidth=2, label='Sheet')
    axes[0, 2].set_xlabel('Time (ns)')
    axes[0, 2].set_ylabel('Content (%)')
    axes[0, 2].set_title('Secondary Structure Formation')
    axes[0, 2].legend()
    axes[0, 2].grid(True, alpha=0.3)
    
    # Native contacts
    axes[1, 0].plot(time_ns, results['native_contacts'], 'm-', linewidth=2)
    axes[1, 0].set_xlabel('Time (ns)')
    axes[1, 0].set_ylabel('Number of Contacts')
    axes[1, 0].set_title('Native Contact Formation')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Torsion cost
    axes[1, 1].plot(time_ns, results['torsion_costs'], 'c-', linewidth=2)
    axes[1, 1].set_xlabel('Time (ns)')
    axes[1, 1].set_ylabel('Average Cost (E_coh)')
    axes[1, 1].set_title('Torsion Angle Optimization')
    axes[1, 1].grid(True, alpha=0.3)
    
    # IR photons
    axes[1, 2].plot(time_ns, results['ir_photons'], 'orange', linewidth=2)
    axes[1, 2].set_xlabel('Time (ns)')
    axes[1, 2].set_ylabel('Cumulative Photons')
    axes[1, 2].set_title('13.8 μm IR Emission')
    axes[1, 2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('mechanism_evolution.png', dpi=150)
    print("\nMechanism evolution plot saved to mechanism_evolution.png")
    plt.close()


def demonstrate_torsion_glyphs(folder: EnhancedThreeLayerFolder):
    """Show the torsion angle glyph distribution"""
    
    # Collect glyph distribution
    glyphs = [r.torsion_state.glyph for r in folder.regions]
    glyph_counts = {i: glyphs.count(i) for i in range(9)}
    
    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Glyph distribution
    ax1.bar(glyph_counts.keys(), glyph_counts.values(), color='purple', alpha=0.7)
    ax1.set_xlabel('Glyph Index')
    ax1.set_ylabel('Count')
    ax1.set_title('Nine-Glyph Distribution')
    ax1.set_xticks(range(9))
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Highlight special glyphs
    ax1.axvspan(3.5, 4.5, alpha=0.2, color='blue', label='Helix (4)')
    ax1.axvspan(-0.5, 0.5, alpha=0.2, color='red', label='Sheet (0)')
    ax1.axvspan(7.5, 8.5, alpha=0.2, color='blue', label='Helix (8)')
    ax1.legend()
    
    # Cost distribution
    costs = [r.torsion_state.cost for r in folder.regions]
    ax2.hist(costs, bins=9, range=(-0.5, 8.5), color='green', alpha=0.7)
    ax2.set_xlabel('Recognition Cost (E_coh)')
    ax2.set_ylabel('Count')
    ax2.set_title('Torsion Cost Distribution')
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('torsion_glyph_analysis.png', dpi=150)
    print("\nTorsion glyph analysis saved to torsion_glyph_analysis.png")
    plt.close()


def main():
    """Run mechanism demonstrations"""
    
    print("\n" + "="*60)
    print("RECOGNITION SCIENCE MECHANISM DEMONSTRATION")
    print("="*60)
    print("\nThis demonstrates our implemented mechanisms without")
    print("waiting for realistic microsecond folding times.")
    
    # Test different sequences
    test_sequences = [
        ("AAAAAAAA", "Poly-alanine (helix-prone)"),
        ("VVVVVVVV", "Poly-valine (sheet-prone)"),
        ("AEAAAKAA", "Helix with breaker"),
    ]
    
    all_results = []
    
    for sequence, description in test_sequences:
        print(f"\n\nTesting: {description}")
        results = demonstrate_mechanisms(sequence, barrier_factor=0.01)  # 100x easier
        all_results.append(results)
        
        # Plot evolution for first sequence
        if len(all_results) == 1:
            plot_mechanism_evolution(results)
            demonstrate_torsion_glyphs(results['folder'])
    
    # Summary
    print("\n" + "="*60)
    print("DEMONSTRATION COMPLETE")
    print("="*60)
    print("\nKey mechanisms demonstrated:")
    print("✓ Information template formation (ps timescale)")
    print("✓ Nine-glyph torsion angle system")
    print("✓ Secondary structure detection")
    print("✓ Native contact formation with phase alignment")
    print("✓ IR photon emission tracking")
    print("✓ Two-timescale physics")
    print("\nNOTE: Barrier crossing was artificially accelerated for demo.")
    print("Real proteins would take microseconds to fold.")


if __name__ == "__main__":
    np.random.seed(42)
    main() 