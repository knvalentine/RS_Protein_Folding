"""
Accelerated Protein Folder with Monte Carlo Barrier Crossing

This module provides efficient simulation of protein folding by:
1. Normal simulation until template formation (ps timescale)
2. Monte Carlo sampling of barrier crossing time
3. Direct jump to folding without simulating waiting period

This allows testing of larger proteins without computational bottleneck.
"""

import numpy as np
from typing import Dict, Optional, Tuple
import time

from enhanced_three_layer_folder import EnhancedThreeLayerFolder

# RS Constants
TAU_0 = 7.33e-15  # s
PHI = (1 + np.sqrt(5)) / 2  # Golden ratio

def calculate_k0_folding(n_residues: int) -> float:
    """
    Calculate size-dependent Arrhenius prefactor.
    
    k₀(n) = (1/8τ₀) × φ^(-n/2) × P_ledger × P_geometric
    
    Args:
        n_residues: Number of residues in protein
        
    Returns:
        k₀ in s^-1
    """
    base_rate = 1 / (8 * TAU_0)  # Eight-beat cycle rate
    phase_factor = PHI ** (-n_residues / 2)  # Phase alignment probability
    ledger_factor = 0.5  # Ledger availability
    geometric_factor = 0.01  # Geometric compatibility
    
    return base_rate * phase_factor * ledger_factor * geometric_factor

class AcceleratedFolder(EnhancedThreeLayerFolder):
    """
    Accelerated protein folder using Monte Carlo for barrier crossing.
    
    Key optimizations:
    1. Monte Carlo sampling of folding initiation time
    2. Option to skip physical folding simulation
    3. Fast template-only mode for information layer studies
    """
    
    def __init__(self, n_residues: int, temperature: float = 310.0,
                 sequence: Optional[str] = None,
                 monte_carlo_folding: bool = True,
                 simulate_physical: bool = True):
        """
        Initialize accelerated folder.
        
        Args:
            n_residues: Number of residues
            temperature: Temperature in K
            sequence: Amino acid sequence
            monte_carlo_folding: Use MC for barrier crossing
            simulate_physical: Whether to simulate physical folding
        """
        super().__init__(n_residues, temperature, sequence)
        
        self.monte_carlo_folding = monte_carlo_folding
        self.simulate_physical = simulate_physical
        self.mc_folding_time = None
        
    def run_until_template(self, max_ps: float = 1000.0) -> Dict[str, float]:
        """
        Run only until template formation (fast).
        
        Args:
            max_ps: Maximum time in picoseconds
            
        Returns:
            Metrics dictionary
        """
        max_ticks = int(max_ps * 1e-12 / TAU_0)
        start_time = time.time()
        
        while self.tick < max_ticks and not self.template_ready:
            self.step()
            
            if self.tick % 10000 == 0:
                current_ps = self.tick * TAU_0 * 1e12
                print(f"  {current_ps:.1f} ps: {len(self.recognition_events)} recognitions")
        
        wall_time = time.time() - start_time
        
        return {
            'template_formed': self.template_ready,
            'template_time_ps': self.template_completion_tick * TAU_0 * 1e12 if self.template_ready else None,
            'recognition_events': len(self.recognition_events),
            'phase_coherence': self.phase_field.get_folding_progress()['coherence'],
            'wall_time_s': wall_time,
            'ticks_simulated': self.tick
        }
    
    def monte_carlo_barrier_crossing(self) -> Tuple[float, int]:
        """
        Use Monte Carlo to sample barrier crossing time.
        
        Returns:
            (crossing_time_us, number_of_attempts)
        """
        # Folding rate at this temperature with size-dependent prefactor
        k0 = calculate_k0_folding(self.n_residues)
        k_fold = k0 * np.exp(-0.18 / self.kT)
        
        # Sample from exponential distribution
        # Time to first success: t = -ln(U) / k
        u = np.random.random()
        crossing_time_s = -np.log(u) / k_fold
        crossing_time_us = crossing_time_s * 1e6
        
        # Estimate number of "attempts" for consistency
        dt = TAU_0  # One attempt per tick
        expected_attempts = int(crossing_time_s / dt)
        
        return crossing_time_us, expected_attempts
    
    def run_accelerated(self, max_us: float = 1000.0,
                       template_timeout_ps: float = 1000.0) -> Dict[str, float]:
        """
        Run accelerated simulation with all optimizations.
        
        Args:
            max_us: Maximum time in microseconds
            template_timeout_ps: Maximum time to wait for template
            
        Returns:
            Complete metrics dictionary
        """
        print(f"\nAccelerated simulation of {self.n_residues} residues")
        print(f"Monte Carlo folding: {self.monte_carlo_folding}")
        print(f"Simulate physical: {self.simulate_physical}")
        print("-" * 50)
        
        start_time = time.time()
        
        # Phase 1: Template formation (normal simulation)
        print("Phase 1: Template formation...")
        template_metrics = self.run_until_template(template_timeout_ps)
        
        if not template_metrics['template_formed']:
            print("✗ Template formation timeout")
            return template_metrics
        
        print(f"✓ Template formed at {template_metrics['template_time_ps']:.1f} ps")
        
        # Phase 2: Barrier crossing
        if self.monte_carlo_folding:
            print("\nPhase 2: Monte Carlo barrier crossing...")
            crossing_time_us, attempts = self.monte_carlo_barrier_crossing()
            self.mc_folding_time = crossing_time_us
            
            # Jump to folding time
            folding_tick = int((crossing_time_us * 1e-6 + 
                               template_metrics['template_time_ps'] * 1e-12) / TAU_0)
            self.tick = folding_tick
            self.folding_initiated = True
            self.barrier_crossing_attempts = attempts
            
            # Calculate k0 and k_fold for display
            k0 = calculate_k0_folding(self.n_residues)
            k_fold = k0 * np.exp(-0.18 / self.kT)
            
            print(f"✓ MC folding at {crossing_time_us:.2f} μs")
            print(f"  k₀ = {k0:.2e} s^-1 (n={self.n_residues})")
            print(f"  k(T) = {k_fold:.2e} s^-1")
        else:
            print("\nPhase 2: Direct barrier crossing simulation...")
            # Original method with acceleration
            max_wait_ticks = int(max_us * 1e-6 / TAU_0)
            acceleration = 1000
            
            while not self.folding_initiated and self.tick < max_wait_ticks:
                self.tick += acceleration
                dt = acceleration * TAU_0
                k0 = calculate_k0_folding(self.n_residues)
                k_fold = k0 * np.exp(-0.18 / self.kT)
                p_initiate = 1 - np.exp(-k_fold * dt)
                
                if np.random.random() < p_initiate:
                    self.folding_initiated = True
                    break
        
        # Phase 3: Physical folding (optional)
        if self.folding_initiated and self.simulate_physical:
            print("\nPhase 3: Physical folding...")
            
            # Simulate a brief period of actual folding
            folding_ticks = min(10000, int(0.1e-6 / TAU_0))  # 0.1 μs or 10k ticks
            
            for i in range(folding_ticks):
                self.step()
                
                if i % 1000 == 0:
                    self._update_torsion_angles()
                    self._detect_secondary_structures()
        
        # Compile final metrics
        final_metrics = self.get_detailed_metrics()
        wall_time = time.time() - start_time
        
        return {
            # Template formation
            'template_formed': True,
            'template_time_ps': template_metrics['template_time_ps'],
            'template_recognitions': template_metrics['recognition_events'],
            
            # Folding initiation
            'folding_initiated': self.folding_initiated,
            'folding_time_us': self.tick * TAU_0 * 1e6,
            'mc_folding_time_us': self.mc_folding_time,
            'barrier_attempts': self.barrier_crossing_attempts,
            
            # Final structure
            'helix_content': final_metrics['helix_content'],
            'sheet_content': final_metrics['sheet_content'],
            'native_contacts': final_metrics['native_contacts'],
            'compactness': final_metrics['compactness'],
            
            # Performance
            'wall_time_s': wall_time,
            'total_ticks': self.tick,
            'acceleration_factor': self.tick / (wall_time / TAU_0) if wall_time > 0 else 0
        }
    
    def estimate_folding_time(self, n_samples: int = 1000) -> Dict[str, float]:
        """
        Estimate folding time distribution using Monte Carlo.
        
        Args:
            n_samples: Number of MC samples
            
        Returns:
            Statistics dictionary
        """
        k0 = calculate_k0_folding(self.n_residues)
        k_fold = k0 * np.exp(-0.18 / self.kT)
        
        # Sample folding times
        samples = -np.log(np.random.random(n_samples)) / k_fold * 1e6  # in μs
        
        return {
            'mean_us': np.mean(samples),
            'median_us': np.median(samples),
            'std_us': np.std(samples),
            'min_us': np.min(samples),
            'max_us': np.max(samples),
            'rate_constant': k_fold,
            'barrier_ev': 0.18,
            'temperature_k': self.temperature
        } 