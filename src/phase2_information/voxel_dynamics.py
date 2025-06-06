"""
Pure Voxel Dynamics for Recognition Science

This module implements discrete voxel transitions instead of
continuous mobility, staying true to RS discrete principles.
"""

import numpy as np
from typing import Tuple, List

# RS Constants
VOXEL_SIZE = 3.35  # Å
PHI = (1 + np.sqrt(5)) / 2
DAMPING_FACTOR = np.sqrt(PHI) / PHI  # For backbone bonds
E_COH = 0.090  # eV
K_B = 8.617e-5  # eV/K

class VoxelDynamics:
    """
    Implements discrete voxel transitions based on recognition events.
    
    Instead of continuous mobility, atoms jump between voxels based on:
    - Local recognition density
    - Information pressure gradient
    - Eight-beat phase alignment
    """
    
    def __init__(self, temperature: float = 310.0):
        self.temperature = temperature
        self.kT = K_B * temperature
        
    def position_to_voxel(self, position: np.ndarray) -> np.ndarray:
        """Convert continuous position to voxel indices"""
        return np.floor(position / VOXEL_SIZE).astype(int)
    
    def voxel_to_position(self, voxel: np.ndarray) -> np.ndarray:
        """Convert voxel indices to center position"""
        return (voxel + 0.5) * VOXEL_SIZE
    
    def compute_transition_probability(self, 
                                     current_voxel: np.ndarray,
                                     target_voxel: np.ndarray,
                                     info_pressure_gradient: np.ndarray,
                                     recognition_density: float,
                                     tick: int) -> float:
        """
        Compute probability of voxel transition based on RS principles.
        
        Args:
            current_voxel: Current voxel indices
            target_voxel: Target voxel indices
            info_pressure_gradient: Information pressure gradient (eV/Å)
            recognition_density: Local recognition event density
            tick: Current simulation tick
            
        Returns:
            Transition probability [0, 1]
        """
        # Check if transition is to adjacent voxel only
        voxel_diff = target_voxel - current_voxel
        if np.abs(voxel_diff).max() > 1:
            return 0.0  # Can only jump to adjacent voxels
        
        # Base probability from recognition density
        # Need sufficient recognitions to enable transition
        base_prob = min(1.0, recognition_density / 1000)  # ~1000 recognitions per transition
        
        # Directional bias from information pressure
        direction = voxel_diff / (np.linalg.norm(voxel_diff) + 1e-10)
        pressure_bias = np.dot(info_pressure_gradient, direction)
        
        # Convert pressure to probability using Boltzmann factor
        if pressure_bias > 0:
            directional_factor = 1.0 - np.exp(-pressure_bias * VOXEL_SIZE / self.kT)
        else:
            directional_factor = 0.0
        
        # Eight-beat modulation
        beat_phase = (tick % 8) * 2 * np.pi / 8
        beat_factor = (1 + np.cos(beat_phase)) / 2
        
        # Backbone damping for connected residues
        damping = DAMPING_FACTOR if np.linalg.norm(voxel_diff) > 0 else 1.0
        
        # Total probability
        prob = base_prob * directional_factor * beat_factor * damping
        
        return min(1.0, prob)
    
    def attempt_voxel_transition(self,
                               position: np.ndarray,
                               info_pressure: np.ndarray,
                               recognition_count: int,
                               tick: int,
                               constrained: bool = False) -> Tuple[np.ndarray, bool]:
        """
        Attempt a voxel transition for an atom.
        
        Args:
            position: Current position (Å)
            info_pressure: Information pressure gradient at position (eV/Å)
            recognition_count: Number of recognition events since last transition
            tick: Current simulation tick
            constrained: Whether atom is constrained (e.g., backbone)
            
        Returns:
            (new_position, transition_occurred)
        """
        current_voxel = self.position_to_voxel(position)
        
        # Consider all 26 adjacent voxels (3D)
        transitions = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dz in [-1, 0, 1]:
                    if dx == 0 and dy == 0 and dz == 0:
                        continue
                    
                    target_voxel = current_voxel + np.array([dx, dy, dz])
                    
                    # Compute transition probability
                    prob = self.compute_transition_probability(
                        current_voxel, target_voxel,
                        info_pressure, recognition_count,
                        tick
                    )
                    
                    if prob > 0:
                        transitions.append((target_voxel, prob))
        
        # Apply backbone constraint
        if constrained:
            # Reduce probabilities for non-preferred directions
            transitions = [(v, p * DAMPING_FACTOR) for v, p in transitions]
        
        # Normalize probabilities
        total_prob = sum(p for _, p in transitions)
        if total_prob == 0:
            return position, False
        
        # Add probability of staying in current voxel
        stay_prob = max(0, 1.0 - total_prob)
        
        # Make stochastic choice
        rand = np.random.random()
        cumulative = 0.0
        
        # Check if we stay
        cumulative += stay_prob
        if rand < cumulative:
            return position, False
        
        # Otherwise, choose transition
        for target_voxel, prob in transitions:
            cumulative += prob / total_prob * (1 - stay_prob)
            if rand < cumulative:
                new_position = self.voxel_to_position(target_voxel)
                return new_position, True
        
        # Shouldn't reach here, but just in case
        return position, False
    
    def compute_voxel_field_energy(self, positions: List[np.ndarray]) -> float:
        """
        Compute total voxel field energy for configuration.
        
        Each occupied voxel contributes E_coh/8 to maintain coherence.
        """
        occupied_voxels = set()
        for pos in positions:
            voxel = tuple(self.position_to_voxel(pos))
            occupied_voxels.add(voxel)
        
        # Energy cost for maintaining voxel coherence
        return len(occupied_voxels) * E_COH / 8 