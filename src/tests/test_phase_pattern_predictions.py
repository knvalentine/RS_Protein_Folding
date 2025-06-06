#!/usr/bin/env python3
"""
Test phase pattern analysis for β-sheet prediction.

This tests whether our new PhasePatternAnalyzer can extract features
from the template that predict final β-sheet content.
"""

import numpy as np
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'core'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from accelerated_folder_v3 import AcceleratedFolderV3
from phase_pattern_analyzer import PhasePatternAnalyzer

# Test proteins with known structures
test_proteins = [
    # (name, sequence, exp_time, temp, expected_beta_registries)
    ('Trp-cage', 'NLYIQWLKDGGPSSGRPPPS', 4.1, 296, 1),  # Has one β-hairpin
    ('BBA5', 'EQYTAKYKGRTFRNEKELRDFIE', 13.0, 298, 1),  # Mixed α/β
    ('Villin', 'LSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF', 0.7, 300, 0),  # All α-helix
    ('WW domain', 'GSKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSG', 13.0, 298, 2),  # β-sheet rich
]

print("Phase Pattern Analysis for β-Sheet Prediction")
print("="*60)

results = []

for name, sequence, exp_time, temp, expected_R in test_proteins:
    print(f"\n{name} ({len(sequence)} residues)")
    print(f"Expected β-registries: {expected_R}")
    print("-"*40)
    
    # Create folder
    folder = AcceleratedFolderV3(
        n_residues=len(sequence),
        sequence=sequence,
        temperature=temp,
        use_template_params=False  # Use original analyzer first
    )
    
    # Run until template forms
    template_result = folder.run_until_template()
    
    # Get phase field
    phase_field = folder.phase_field
    
    # Original analysis
    if hasattr(folder, 'template_analysis') and folder.template_analysis:
        original_params = folder.get_parameter_summary()
    else:
        # Manually analyze if needed
        from pattern_analyzer import PatternAnalyzer
        analyzer = PatternAnalyzer()
        positions = folder.positions.copy()
        torsion_states = getattr(folder, 'torsion_states', 
                                np.random.randint(0, 9, folder.n_residues))
        template_analysis = analyzer.analyze_template(
            phase_field, positions, torsion_states
        )
        # Extract beta registries from the analysis
        beta_registries = folder.pattern_analyzer._count_beta_registries(
            None, {}, torsion_states, positions
        )
        original_params = {
            'beta_registries': beta_registries,
            'P_geom': template_analysis.p_geom,
            'loops': template_analysis.n_loops,
            'contact_order': template_analysis.contact_order,
            'helix_axes': 0  # Not directly stored
        }
    original_R = original_params.get('beta_registries', 0)
    
    # New phase pattern analysis
    analyzer = PhasePatternAnalyzer(phase_field)
    features = analyzer.get_predictive_features()
    predicted_R = features['predicted_beta_registries']
    
    print(f"Original R detection: {original_R}")
    print(f"Phase pattern prediction: {predicted_R}")
    print(f"Actual R: {expected_R}")
    
    # Show key metrics
    print(f"\nPhase metrics:")
    print(f"  Long-range correlation: {features['long_range_correlation']:.3f}")
    print(f"  Phase frustration: {features['phase_frustration']:.3f}")
    print(f"  Channel coupling: {features['channel_coupling']:.3f}")
    print(f"  Max coherence: {features['max_coherence']:.3f}")
    
    # Calculate folding time with new prediction
    if predicted_R != original_R:
        # Recalculate P_geom with new R value
        import math
        PHI = (1 + math.sqrt(5)) / 2
        
        # Get other parameters
        loops = original_params.get('loops', 0)
        CO = original_params.get('contact_order', 0)
        L = original_params.get('helix_axes', 0)
        
        # New P_geom with predicted R
        new_P_geom = (PHI**(-loops/2) * 
                      PHI**(-math.sqrt(PHI)*CO) * 
                      PHI**(L/2) * 
                      PHI**(-predicted_R))
        
        # Cap at reasonable value
        new_P_geom = min(new_P_geom, 0.1)
        
        print(f"\nOriginal P_geom: {original_params.get('P_geom', 1.0):.6f}")
        print(f"New P_geom with predicted R: {new_P_geom:.6f}")
        
        # Update parameters and run
        folder.use_template_params = True
        folder.template_params = original_params.copy()
        folder.template_params['beta_registries'] = predicted_R
        folder.template_params['P_geom'] = new_P_geom
        
    # Run folding
    result = folder.run_accelerated(max_us=10000.0)
    
    if result.get('total_time_us'):
        total_time = result['total_time_us']
        ratio = total_time / exp_time
        
        print(f"\nFolding results:")
        print(f"RS prediction: {total_time:.1f} μs")
        print(f"Experimental: {exp_time} μs")
        print(f"Ratio: {ratio:.2f}")
        
        if 0.2 <= ratio <= 5.0:
            print("✅ Within 5× of experiment!")
        else:
            print(f"❌ Off by {max(ratio, 1/ratio):.0f}×")
        
        results.append({
            'protein': name,
            'predicted_R': predicted_R,
            'actual_R': expected_R,
            'ratio': ratio
        })

print("\n" + "="*60)
print("Summary:")
print("-"*60)

# Check if predictions improved
correct_predictions = sum(1 for r in results if r['predicted_R'] == r['actual_R'])
total = len(results)
print(f"β-registry predictions: {correct_predictions}/{total} correct")

# Check folding time accuracy
within_5x = sum(1 for r in results if 0.2 <= r['ratio'] <= 5.0)
print(f"Folding times within 5×: {within_5x}/{total}")

# Average error
avg_error = np.mean([max(r['ratio'], 1/r['ratio']) for r in results])
print(f"Average fold error: {avg_error:.1f}×") 