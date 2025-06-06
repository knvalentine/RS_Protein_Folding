"""
Pattern Analyzer V2 - Reading the Full Phase Information

This enhanced analyzer extracts not just topology but the actual phase
relationships that determine folding kinetics. Based on our deeper
understanding that folding is information organizing matter.
"""

import numpy as np
from typing import Dict, List, Tuple, Set, Optional
from dataclasses import dataclass
from collections import defaultdict
import networkx as nx

# RS constants
PHI = 1.618033988749895
E_COH = 0.090  # eV
TAU_0 = 7.33e-15  # s

@dataclass
class EnhancedTemplateAnalysis:
    """Results from deep phase analysis of template."""
    # Original topological features
    barrier_coins: int
    barrier_ev: float
    p_ledger: float
    p_geom: float
    n_voxels: int
    n_components: int
    n_loops: int
    
    # New phase-based features
    phase_coherence_length: float  # ξ in nm
    information_flow_rate: float   # bits/ps
    recognition_density: float     # events/voxel
    phase_frustration: float       # 0-1 scale
    channel_occupancy: List[float] # usage of 8 channels
    effective_barrier_ev: float    # barrier from phase analysis
    
    
class PatternAnalyzerV2:
    """Enhanced analyzer that reads phase information, not just topology."""
    
    def __init__(self, voxel_size: float = 0.335e-9):
        self.voxel_size = voxel_size
        
    def analyze_template(self, phase_field, positions: np.ndarray, 
                        torsion_states: Optional[np.ndarray] = None) -> EnhancedTemplateAnalysis:
        """
        Extract both topological and phase information from template.
        
        Args:
            phase_field: PhasePatternField after template formation
            positions: Current atom positions [n_atoms, 3]
            torsion_states: Optional torsion states (will infer if not provided)
            
        Returns:
            EnhancedTemplateAnalysis with full phase information
        """
        # Build voxel graph
        voxel_graph, voxel_map = self._build_voxel_graph(phase_field, positions)
        
        # Original topological analysis
        n_components = nx.number_connected_components(voxel_graph)
        n_loops = self._count_independent_loops(voxel_graph)
        n_voxels = voxel_graph.number_of_nodes()
        
        # NEW: Extract phase information
        phase_coherence_length = self._calculate_coherence_length(phase_field)
        info_flow_rate = self._calculate_information_flow(phase_field)
        recognition_density = self._calculate_recognition_density(phase_field)
        phase_frustration = self._calculate_phase_frustration(phase_field, voxel_graph)
        channel_occupancy = self._analyze_channel_usage(phase_field)
        
        # NEW: Derive barrier from phase properties, not categories
        effective_barrier = self._calculate_effective_barrier(
            phase_coherence_length, info_flow_rate, 
            recognition_density, phase_frustration
        )
        
        # NEW: Derive geometric factors from phase flow
        p_ledger = self._calculate_p_ledger_from_channels(channel_occupancy)
        p_geom = self._calculate_p_geom_from_coherence(
            phase_coherence_length, phase_frustration
        )
        
        # For compatibility, still calculate integer coins
        barrier_coins = int(np.round(effective_barrier / E_COH))
        barrier_coins = max(2, min(4, barrier_coins))  # Keep in reasonable range
        
        return EnhancedTemplateAnalysis(
            barrier_coins=barrier_coins,
            barrier_ev=barrier_coins * E_COH,
            p_ledger=p_ledger,
            p_geom=p_geom,
            n_voxels=n_voxels,
            n_components=n_components,
            n_loops=n_loops,
            phase_coherence_length=phase_coherence_length,
            information_flow_rate=info_flow_rate,
            recognition_density=recognition_density,
            phase_frustration=phase_frustration,
            channel_occupancy=channel_occupancy,
            effective_barrier_ev=effective_barrier
        )
    
    def _build_voxel_graph(self, phase_field, positions: np.ndarray) -> Tuple[nx.Graph, Dict]:
        """Build voxel connectivity graph (same as before)."""
        voxel_indices = phase_field.get_voxel_indices(positions)
        
        voxel_map = defaultdict(list)
        for i, vox_idx in enumerate(voxel_indices):
            if vox_idx is not None:
                voxel_map[vox_idx].append(i)
        
        G = nx.Graph()
        occupied_voxels = list(voxel_map.keys())
        G.add_nodes_from(occupied_voxels)
        
        for i, vox1 in enumerate(occupied_voxels):
            for vox2 in occupied_voxels[i+1:]:
                if self._are_adjacent_voxels(vox1, vox2, phase_field.grid_size):
                    G.add_edge(vox1, vox2)
        
        return G, dict(voxel_map)
    
    def _are_adjacent_voxels(self, vox1: Tuple, vox2: Tuple, grid_size: np.ndarray) -> bool:
        """Check if two voxel indices are adjacent."""
        diff = np.abs(np.array(vox1) - np.array(vox2))
        return np.all(diff <= 1) and np.sum(diff) > 0
    
    def _count_independent_loops(self, graph: nx.Graph) -> int:
        """Count independent loops in graph."""
        if graph.number_of_edges() == 0:
            return 0
        
        n_loops = 0
        for component in nx.connected_components(graph):
            subgraph = graph.subgraph(component)
            n_loops += subgraph.number_of_edges() - subgraph.number_of_nodes() + 1
        
        return n_loops
    
    def _calculate_coherence_length(self, phase_field) -> float:
        """
        Calculate phase coherence length ξ.
        This tells us how far phase information propagates.
        """
        # Get phase pattern data
        phase_pattern = phase_field.phase_pattern  # [n_residues, 8_channels]
        n_residues = phase_pattern.shape[0]
        
        if n_residues < 2:
            return self.voxel_size * 1e9  # nm
        
        # Calculate phase-phase correlation vs distance
        correlations = []
        distances = []
        
        for i in range(n_residues):
            for j in range(i+1, min(i+20, n_residues)):  # Look up to 20 residues away
                # Phase correlation across all channels
                phi_i = phase_pattern[i, :]
                phi_j = phase_pattern[j, :]
                
                # Complex correlation
                correlation = np.abs(np.mean(np.exp(1j * (phi_i - phi_j))))
                correlations.append(correlation)
                distances.append((j - i) * self.voxel_size * 1e9)  # nm
        
        if not correlations:
            return self.voxel_size * 1e9
        
        # Fit exponential decay to get coherence length
        # C(r) = exp(-r/ξ)
        correlations = np.array(correlations)
        distances = np.array(distances)
        
        # Simple estimate: distance where correlation drops to 1/e
        threshold_idx = np.where(correlations < 1/np.e)[0]
        if len(threshold_idx) > 0:
            xi = distances[threshold_idx[0]]
        else:
            xi = distances[-1]  # Coherence extends beyond our measurement
        
        return xi
    
    def _calculate_information_flow(self, phase_field) -> float:
        """
        Calculate information flow rate through the template.
        High flow = fast folding.
        """
        phase_pattern = phase_field.phase_pattern
        
        # Information is phase gradient magnitude
        info_gradients = []
        
        for channel in range(8):
            phases = phase_pattern[:, channel]
            if len(phases) > 1:
                # Discrete gradient
                grad = np.diff(phases)
                info_gradients.extend(np.abs(grad))
        
        if not info_gradients:
            return 0.0
        
        # Average information flow in radians per voxel
        avg_gradient = np.mean(info_gradients)
        
        # Convert to bits/ps using recognition frequency
        # Information = gradient * frequency * log2(8 channels)
        info_rate = avg_gradient * (1/(8*TAU_0)) * 1e-12 * 3  # bits/ps
        
        return info_rate
    
    def _calculate_recognition_density(self, phase_field) -> float:
        """
        Calculate density of recognition events.
        Dense recognition = complex folding.
        """
        # Recognition events are stored in the recognition matrix
        recognition_matrix = phase_field.recognition_matrix
        
        # Count total recognition events
        n_events = np.sum(recognition_matrix) / 2  # Matrix is symmetric
        
        # Normalize by number of voxels
        n_residues = recognition_matrix.shape[0]
        density = n_events / n_residues if n_residues > 0 else 0
        
        return density
    
    def _calculate_phase_frustration(self, phase_field, voxel_graph) -> float:
        """
        Calculate phase frustration from triangular loops.
        Frustrated phases can't all be satisfied simultaneously.
        """
        frustration = 0.0
        n_triangles = 0
        
        # Find all triangles in the voxel graph
        for node in voxel_graph.nodes():
            neighbors = list(voxel_graph.neighbors(node))
            for i in range(len(neighbors)):
                for j in range(i+1, len(neighbors)):
                    if voxel_graph.has_edge(neighbors[i], neighbors[j]):
                        # Found a triangle
                        n_triangles += 1
                        
                        # Get phases from voxel edges
                        edges = phase_field.voxel_edges
                        
                        # Map voxel indices to residue indices (simplified)
                        # In reality, would need proper mapping
                        phase_sum = 0
                        for edge in [(node, neighbors[i]), 
                                   (neighbors[i], neighbors[j]),
                                   (neighbors[j], node)]:
                            if edge in edges:
                                phase_sum += edges[edge]
                            elif (edge[1], edge[0]) in edges:
                                phase_sum -= edges[(edge[1], edge[0])]
                        
                        # Frustration is deviation from 0 or 2π
                        frustration += min(abs(phase_sum), abs(phase_sum - 2*np.pi))
        
        # Normalize to 0-1 scale
        if n_triangles > 0:
            frustration = frustration / (n_triangles * np.pi)
            frustration = min(1.0, frustration)
        
        return frustration
    
    def _analyze_channel_usage(self, phase_field) -> List[float]:
        """
        Analyze how much each of the 8 channels is used.
        Proteins using all channels have different folding.
        """
        phase_pattern = phase_field.phase_pattern
        channel_usage = []
        
        for channel in range(8):
            # RMS amplitude in this channel
            channel_phases = phase_pattern[:, channel]
            rms = np.sqrt(np.mean(channel_phases**2))
            channel_usage.append(rms)
        
        # Normalize so max channel = 1.0
        max_usage = max(channel_usage) if channel_usage else 1.0
        if max_usage > 0:
            channel_usage = [u/max_usage for u in channel_usage]
        
        return channel_usage
    
    def _calculate_effective_barrier(self, coherence_length: float,
                                   info_flow: float, 
                                   rec_density: float,
                                   frustration: float) -> float:
        """
        Return the universal Recognition-Science barrier ΔC = 2 E_COH = 0.18 eV.
        RS principles forbid protein-specific barrier adjustments; topology and
        phase only tune the Arrhenius pre-factor (k₀), never ΔC itself.
        """
        return 2 * E_COH
    
    def _calculate_p_ledger_from_channels(self, channel_occupancy: List[float]) -> float:
        """
        Ledger availability from channel usage pattern.
        Using all 8 channels = more ledger coordination needed.
        """
        # Count significantly used channels (> 0.2 of max)
        active_channels = sum(1 for c in channel_occupancy if c > 0.2)
        
        # Each active channel beyond the first reduces availability
        p_ledger = PHI ** (-(active_channels - 1) / 4)
        
        return p_ledger
    
    def _calculate_p_geom_from_coherence(self, coherence_length: float,
                                       frustration: float) -> float:
        """
        Geometric factor from phase coherence properties.
        Long coherence + low frustration = favorable geometry.
        """
        # Base geometric factor
        p_geom = 0.01
        
        # Coherence bonus/penalty
        # Long coherence (>10 nm) is favorable
        if coherence_length > 10:
            p_geom *= PHI ** ((coherence_length - 10) / 20)
        elif coherence_length < 5:
            p_geom *= PHI ** (-(5 - coherence_length) / 5)
        
        # Frustration penalty
        p_geom *= (1 - frustration)
        
        # Keep in reasonable bounds
        p_geom = np.clip(p_geom, 1e-4, 1.0)
        
        return p_geom 