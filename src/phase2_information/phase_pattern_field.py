"""
Phase Pattern Field - Information Layer for Recognition Science

This module implements the information layer that bridges quantum recognition
events and physical protein folding. Phase patterns form in ~65 picoseconds
and guide microsecond-scale physical reconfiguration.

Key principles:
- 8-channel phase architecture (from RS eight-beat cycle)
- No forces or empirical parameters
- Information precedes physical change
- Golden ratio (φ) scaling throughout
"""

import numpy as np
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass

# Recognition Science constants
TAU_0 = 7.33e-15  # s (fundamental tick)
E_COH = 0.090  # eV (one recognition quantum)
PHI = (1 + np.sqrt(5)) / 2  # Golden ratio
RECOGNITION_DISTANCE = 6.5  # Å

# Eight-channel architecture
N_CHANNELS = 8
CHANNEL_PHASES = np.array([i * 137.5 for i in range(N_CHANNELS)]) * np.pi / 180  # Golden angle spacing

@dataclass
class RecognitionEvent:
    """A discrete recognition event between two residues"""
    residue_i: int
    residue_j: int
    phase_shift: float  # radians
    tick: int
    energy: float = E_COH  # Always one quantum


class PhasePatternField:
    """
    Information field that accumulates phase patterns from recognition events.
    
    The field stores phase information for each residue across 8 channels,
    representing the information template that will guide physical folding.
    """
    
    def __init__(self, n_residues: int):
        """
        Initialize phase pattern field for a protein.
        
        Args:
            n_residues: Number of amino acid residues
        """
        self.n_residues = n_residues
        
        # 8-channel phase storage for each residue
        # Shape: (n_residues, 8)
        self.phase_pattern = np.zeros((n_residues, N_CHANNELS))
        
        # Track which residues have been recognized
        self.recognition_matrix = np.zeros((n_residues, n_residues), dtype=bool)
        
        # Voxel graph encoding folding pathway
        self.voxel_edges: Dict[Tuple[int, int], float] = {}
        
        # Track total recognition events
        self.recognition_count = 0
        self.tick_count = 0
        
        # Coherence metrics
        self.coherence_time = None
        self.is_template_complete = False
        
    def add_recognition(self, event: RecognitionEvent) -> None:
        """
        Process a recognition event and update phase patterns.
        
        Recognition events transfer phase information between residues,
        building up the 8-channel pattern that encodes the folding pathway.
        
        Args:
            event: RecognitionEvent containing residue pair and phase info
        """
        i, j = event.residue_i, event.residue_j
        
        # Mark recognition in matrix
        self.recognition_matrix[i, j] = True
        self.recognition_matrix[j, i] = True
        
        # Update voxel graph
        self.voxel_edges[(i, j)] = event.phase_shift
        
        # Distribute phase shift across channels
        # Each channel gets phase modulated by its characteristic angle
        for channel in range(N_CHANNELS):
            channel_phase = CHANNEL_PHASES[channel]
            
            # Phase coupling strength follows golden ratio
            coupling = np.cos(event.phase_shift - channel_phase) / PHI
            
            # Update phase pattern (conserving total phase)
            self.phase_pattern[i, channel] += coupling
            self.phase_pattern[j, channel] -= coupling
        
        self.recognition_count += 1
        self.tick_count = event.tick
        
        # Check for template completion
        if not self.is_template_complete:
            self._check_coherence()
    
    def _check_coherence(self) -> None:
        """
        Check if phase pattern has achieved coherence.
        
        Coherence requires:
        1. Phase balance across channels (ledger neutrality)
        2. Spanning tree in recognition network
        3. Sufficient phase amplitude
        """
        # Check 1: Phase balance (ledger must balance)
        channel_sums = np.sum(self.phase_pattern, axis=0)
        phase_balanced = np.allclose(channel_sums, 0, atol=1e-10)
        
        # Check 2: Network connectivity (information must span protein)
        connected = self._is_connected()
        
        # Check 3: Phase amplitude (information must be strong enough)
        max_amplitude = np.max(np.abs(self.phase_pattern))
        sufficient_amplitude = max_amplitude > 1.0  # In units of E_COH
        
        if phase_balanced and connected and sufficient_amplitude:
            self.is_template_complete = True
            self.coherence_time = self.tick_count * TAU_0
            
    def _is_connected(self) -> bool:
        """Check if recognition network forms a spanning tree"""
        if self.n_residues <= 1:
            return True
            
        # Simple DFS to check connectivity
        visited = np.zeros(self.n_residues, dtype=bool)
        stack = [0]
        visited[0] = True
        count = 1
        
        while stack and count < self.n_residues:
            current = stack.pop()
            for neighbor in range(self.n_residues):
                if self.recognition_matrix[current, neighbor] and not visited[neighbor]:
                    visited[neighbor] = True
                    stack.append(neighbor)
                    count += 1
                    
        return count == self.n_residues
    
    def compute_information_pressure(self, positions: np.ndarray) -> np.ndarray:
        """
        Compute information pressure gradient for physical update.
        
        The mismatch between the information template (phase pattern) and
        current physical configuration creates pressure for reconfiguration.
        
        Args:
            positions: Current 3D positions of residues, shape (n_residues, 3)
            
        Returns:
            pressure_gradient: Information pressure on each residue, shape (n_residues, 3)
        """
        pressure_gradient = np.zeros_like(positions)
        
        # For each recognition edge in the voxel graph
        for (i, j), target_phase in self.voxel_edges.items():
            # Current physical separation
            r_ij = positions[j] - positions[i]
            distance = np.linalg.norm(r_ij)
            
            if distance < 1e-10:  # Avoid division by zero
                continue
                
            # Information wants specific phase relationship
            # Map phase to preferred distance using golden ratio
            target_distance = RECOGNITION_DISTANCE * (PHI ** (target_phase / (2 * np.pi)))
            
            # Pressure proportional to distance mismatch
            pressure_magnitude = (distance - target_distance) * E_COH
            
            # Apply pressure along connection vector
            pressure_direction = r_ij / distance
            pressure_gradient[i] += pressure_magnitude * pressure_direction
            pressure_gradient[j] -= pressure_magnitude * pressure_direction
            
        return pressure_gradient
    
    def get_phase_coherence_metric(self) -> float:
        """
        Calculate overall phase coherence (0 to 1).
        
        Returns:
            coherence: 1.0 = perfect coherence, 0.0 = no coherence
        """
        if self.recognition_count == 0:
            return 0.0
            
        # Coherence based on phase alignment across channels
        coherence = 0.0
        
        for channel in range(N_CHANNELS):
            channel_phases = self.phase_pattern[:, channel]
            
            # Calculate phase order parameter (like Kuramoto model)
            if np.any(channel_phases):
                z = np.mean(np.exp(1j * channel_phases))
                coherence += np.abs(z)
                
        return coherence / N_CHANNELS
    
    def get_folding_progress(self) -> Dict[str, float]:
        """
        Get metrics describing folding progress.
        
        Returns:
            Dictionary with progress metrics
        """
        return {
            'recognition_count': self.recognition_count,
            'time_ps': self.tick_count * TAU_0 * 1e12,
            'coherence': self.get_phase_coherence_metric(),
            'template_complete': self.is_template_complete,
            'coherence_time_ps': self.coherence_time * 1e12 if self.coherence_time else None,
            'network_density': np.sum(self.recognition_matrix) / (self.n_residues * (self.n_residues - 1)),
        }
    
    def get_channel_summary(self) -> np.ndarray:
        """
        Get summary of phase amplitude in each channel.
        
        Returns:
            Array of RMS phase amplitudes for each channel
        """
        rms_amplitudes = np.zeros(N_CHANNELS)
        for channel in range(N_CHANNELS):
            rms_amplitudes[channel] = np.sqrt(np.mean(self.phase_pattern[:, channel]**2))
        return rms_amplitudes
    
    def get_voxel_indices(self, positions: np.ndarray) -> List[Optional[Tuple[int, int, int]]]:
        """
        Map 3D positions to voxel indices.
        
        Args:
            positions: Array of 3D positions, shape (n_atoms, 3)
            
        Returns:
            List of voxel indices (i, j, k) or None if outside grid
        """
        # Define voxel grid parameters
        voxel_size = 3.35  # Å (0.335 nm)
        
        # Find bounding box
        min_pos = np.min(positions, axis=0)
        max_pos = np.max(positions, axis=0)
        
        # Add padding
        padding = 2 * voxel_size
        min_pos -= padding
        max_pos += padding
        
        # Calculate grid size
        grid_size = np.ceil((max_pos - min_pos) / voxel_size).astype(int)
        self.grid_size = grid_size  # Store for pattern analyzer
        
        # Map positions to voxel indices
        voxel_indices = []
        for pos in positions:
            # Calculate voxel index
            voxel_idx = np.floor((pos - min_pos) / voxel_size).astype(int)
            
            # Check if within grid
            if np.all(voxel_idx >= 0) and np.all(voxel_idx < grid_size):
                voxel_indices.append(tuple(voxel_idx))
            else:
                voxel_indices.append(None)
                
        return voxel_indices 