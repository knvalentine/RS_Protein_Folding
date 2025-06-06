"""
Accelerated Folder V3 with Template-Driven Parameters

This version extracts parameters from the formed template rather than
guessing from sequence. This is the RS-pure approach: let the emerged
pattern tell us what it is, don't try to predict it.
"""

import numpy as np
from typing import Dict, Optional, Tuple
import time

from accelerated_folder import AcceleratedFolder
from pattern_analyzer import PatternAnalyzer, TemplateAnalysis
from enhanced_three_layer_folder import EnhancedThreeLayerFolder

# Golden ratio
PHI = (1 + np.sqrt(5)) / 2
TAU_0 = 7.33e-15  # s


class AcceleratedFolderV3(AcceleratedFolder):
    """
    Template-driven accelerated folder.
    
    Key innovation: Run until template forms, analyze the pattern,
    then use extracted parameters for Monte Carlo folding.
    """
    
    def __init__(self, n_residues: int, temperature: float = 310.0,
                 sequence: Optional[str] = None,
                 monte_carlo_folding: bool = True,
                 simulate_physical: bool = True,
                 use_template_params: bool = True):
        """
        Initialize with option for template-driven parameters.
        
        Args:
            use_template_params: Whether to extract params from template
        """
        super().__init__(n_residues, temperature, sequence,
                        monte_carlo_folding, simulate_physical)
        
        self.use_template_params = use_template_params
        self.pattern_analyzer = PatternAnalyzer()
        self.template_analysis = None
        
    def run_accelerated(self, max_us: float = 1000.0) -> Dict:
        """
        Run with template analysis if enabled.
        """
        if not self.use_template_params:
            # Use parent implementation
            return super().run_accelerated(max_us)
        
        # Phase 1: Run until template forms
        print("Phase 1: Template formation...")
        template_metrics = self._run_until_template_with_analysis()
        
        if not template_metrics['template_formed']:
            print("Warning: Template did not form, using default parameters")
            return super().run_accelerated(max_us)
        
        # Phase 2: Monte Carlo with template-derived parameters
        print("\nPhase 2: Monte Carlo folding with template parameters...")
        
        if self.monte_carlo_folding:
            mc_time_us, mc_attempts = self.monte_carlo_barrier_crossing()
            total_time_us = template_metrics['template_time_us'] + mc_time_us
            
            metrics = {
                'template_time_us': template_metrics['template_time_us'],
                'mc_folding_time_us': mc_time_us,
                'total_time_us': total_time_us,
                'mc_attempts': mc_attempts,
                'template_analysis': self.template_analysis,
                'parameters_used': self.get_parameter_summary()
            }
            
            print(f"\nTemplate formed in: {template_metrics['template_time_us']:.3f} μs")
            print(f"MC folding time: {mc_time_us:.1f} μs")
            print(f"Total time: {total_time_us:.1f} μs")
            
            return metrics
        else:
            # Just return template metrics
            return template_metrics
    
    def _run_until_template_with_analysis(self) -> Dict:
        """
        Run until template forms and analyze the pattern.
        """
        # Use parent class's run_until_template method
        template_metrics = super().run_until_template()
        
        if not template_metrics['template_formed']:
            return {
                'template_formed': False,
                'template_time_us': 0.0
            }
        
        # Get the final state from the folder
        positions = self.positions.copy()
        # Initialize torsion states randomly if not tracked
        # This avoids the all-zeros = all-sheet problem
        torsion_states = getattr(self, 'torsion_states', 
                                np.random.randint(0, 9, self.n_residues))
        phase_field = self.phase_field
        
        # Analyze the template
        self.template_analysis = self.pattern_analyzer.analyze_template(
            phase_field, positions, torsion_states
        )
        
        # Store the extracted parameters
        self._apply_template_parameters(self.template_analysis)
        
        return {
            'template_formed': True,
            'template_time_us': template_metrics['template_time_ps'] / 1000.0,  # ps to μs
            'template_analysis': self.template_analysis
        }
    
    def _apply_template_parameters(self, analysis: TemplateAnalysis):
        """Apply parameters extracted from template."""
        self.barrier_ev = analysis.barrier_ev
        self.barrier_coins = analysis.barrier_coins
        self.p_ledger = analysis.p_ledger
        self.p_geom = analysis.p_geom
        self.path_entropy = analysis.path_entropy
        self.mobility_anisotropy = analysis.mobility_anisotropy
        
        # Store additional info
        self.n_voxels = analysis.n_voxels
        self.n_components = analysis.n_components
        self.n_loops = analysis.n_loops
        self.contact_order = analysis.contact_order
        
    def monte_carlo_barrier_crossing(self) -> Tuple[float, int]:
        """
        Monte Carlo with template-derived parameters.
        """
        if self.use_template_params and self.template_analysis:
            # Use template-derived parameters
            barrier = self.barrier_ev
            
            # Calculate k0 with template-derived factors
            base_k0 = 1 / (8 * TAU_0) * PHI ** (-self.n_residues / 2)
            # Path entropy: 0 = unique paths (good), 1 = redundant (bad)
            # So we use (1 - path_entropy) as the efficiency factor
            path_efficiency = 1.0 - self.path_entropy
            # Mobility: higher anisotropy reduces efficiency
            mobility_factor = 1.0 / (1.0 + self.mobility_anisotropy)
            k0 = base_k0 * self.p_ledger * self.p_geom * path_efficiency * mobility_factor
        else:
            # Fall back to defaults
            barrier = 0.18
            from accelerated_folder import calculate_k0_folding
            k0 = calculate_k0_folding(self.n_residues)
        
        # Folding rate
        k_fold = k0 * np.exp(-barrier / self.kT)
        
        # Sample from exponential distribution
        u = np.random.random()
        crossing_time_s = -np.log(u) / k_fold
        crossing_time_us = crossing_time_s * 1e6
        
        # Estimate attempts
        dt = TAU_0
        expected_attempts = int(crossing_time_s / dt)
        
        # Store for analysis
        self.k0_used = k0
        self.barrier_used = barrier
        self.k_fold = k_fold
        
        return crossing_time_us, expected_attempts
    
    def get_parameter_summary(self) -> Dict:
        """Get summary of parameters used."""
        
        if not self.use_template_params or not self.template_analysis:
            return {
                'mode': 'default',
                'barrier_ev': 0.18,
                'p_ledger': 0.5,
                'p_geom': 0.01
            }
        
        analysis = self.template_analysis
        return {
            'mode': 'template-driven',
            'barrier_ev': analysis.barrier_ev,
            'barrier_coins': analysis.barrier_coins,
            'p_ledger': analysis.p_ledger,
            'p_geom': analysis.p_geom,
            'path_entropy': analysis.path_entropy,
            'mobility_anisotropy': analysis.mobility_anisotropy,
            'n_voxels': analysis.n_voxels,
            'n_components': analysis.n_components,
            'n_loops': analysis.n_loops,
            'unique_rungs': len(analysis.unique_rungs),
            'contact_order': analysis.contact_order,
            'helix_fraction': analysis.helix_fraction,
            'sheet_fraction': analysis.sheet_fraction
        }


def test_template_driven_folding():
    """Test the template-driven approach on our problem proteins."""
    
    test_proteins = [
        ('Trp-cage', 'NLYIQWLKDGGPSSGRPPPS', 4.1, 296),
        ('Villin', 'LSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF', 0.7, 300),
        ('BBA5', 'EQYTAKYKGRTFRNEKELRDFIE', 13.0, 298),
        ('WW domain', 'GSKLPPGWEKRMSRSSGRVYYFNHITNASQWERPS', 13.0, 298)
    ]
    
    print("="*80)
    print("TEMPLATE-DRIVEN PARAMETER EXTRACTION TEST")
    print("="*80)
    
    for name, sequence, exp_time_us, exp_temp_k in test_proteins:
        print(f"\n{'='*60}")
        print(f"{name} - Template Analysis")
        print(f"{'='*60}")
        print(f"Sequence length: {len(sequence)}")
        print(f"Experimental: {exp_time_us} μs at {exp_temp_k}K")
        
        # Create folder
        folder = AcceleratedFolderV3(
            n_residues=len(sequence),
            temperature=exp_temp_k,
            sequence=sequence,
            use_template_params=True
        )
        
        # Run with template analysis
        metrics = folder.run_accelerated(max_us=10000.0)
        
        if metrics.get('template_analysis'):
            print("\nTemplate Analysis Results:")
            params = metrics['parameters_used']
            print(f"  Voxels: {params['n_voxels']}")
            print(f"  Components: {params['n_components']}")
            print(f"  Loops: {params['n_loops']}")
            print(f"  Unique rungs: {params['unique_rungs']}")
            print(f"  Contact order: {params['contact_order']:.3f}")
            print(f"  Helix fraction: {params['helix_fraction']:.1%}")
            print(f"  Sheet fraction: {params['sheet_fraction']:.1%}")
            
            print(f"\nExtracted Parameters:")
            print(f"  Barrier: {params['barrier_coins']} coins ({params['barrier_ev']:.3f} eV)")
            print(f"  P_ledger: {params['p_ledger']:.3f}")
            print(f"  P_geom: {params['p_geom']:.4f}")
            print(f"  Path entropy: {params['path_entropy']:.4f}")
            print(f"  Mobility anisotropy: {params['mobility_anisotropy']:.4f}")
            
            print(f"\nFolding Results:")
            print(f"  Template time: {metrics['template_time_us']:.3f} μs")
            print(f"  MC folding time: {metrics['mc_folding_time_us']:.1f} μs")
            print(f"  Total RS time: {metrics['total_time_us']:.1f} μs")
            print(f"  Ratio (RS/Exp): {metrics['total_time_us']/exp_time_us:.2f}")
            
            # Check if we're within order of magnitude
            ratio = metrics['total_time_us'] / exp_time_us
            if 0.1 <= ratio <= 10:
                print("  ✅ Within order of magnitude!")
            else:
                print(f"  ❌ Off by {max(ratio, 1/ratio):.1f}x")
    
    print("\n" + "="*80)
    print("SUMMARY: Template-driven parameters extracted from emerged patterns")
    print("This is the RS way - let reality tell us what it is!")
    print("="*80)


if __name__ == "__main__":
    np.random.seed(42)
    test_template_driven_folding() 