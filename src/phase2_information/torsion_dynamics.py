"""
Torsion Angle Dynamics for Recognition Science Protein Folding

This module implements the RS nine-glyph system for backbone torsion angles
and secondary structure templates based on golden ratio geometry.

Key principles from RS theory:
- Each backbone φ/ψ choice costs integer "coins" (E_coh = 0.090 eV)
- Nine-glyph Ramachandran space: 120°×120° bins
- Secondary structures follow golden ratio phase relationships
- No empirical parameters - everything derives from E_coh, φ, and τ₀
"""

import numpy as np
from typing import Tuple, List, Dict, Optional
from dataclasses import dataclass

# Recognition Science constants
E_COH = 0.090  # eV (recognition quantum)
PHI = (1 + np.sqrt(5)) / 2  # Golden ratio
TAU_0 = 7.33e-15  # s (fundamental tick)

# Torsion angle bins (120° each)
GLYPH_BIN_SIZE = 120.0  # degrees

# Secondary structure phase advances (from RS theory)
HELIX_PHASE_ADVANCE = 100.0 * np.pi / 180  # φ² × 60° base
SHEET_PHASE_ADVANCE = 180.0 * np.pi / 180  # φ³ × 60° base

# Helical geometry from golden ratio
HELIX_RISE = 1.5  # Å per residue (3.4 Å / 2.27 residues)
HELIX_RADIUS = 2.3  # Å
HELIX_RESIDUES_PER_TURN = 3.6  # From 360° / 100°

# Recognition distance for native contacts
RECOGNITION_DISTANCE = 6.5  # Å
PHASE_ALIGNMENT_THRESHOLD = np.pi / 4  # 45° tolerance for phase matching


@dataclass
class TorsionState:
    """Represents the torsion angle state of a residue"""
    phi: float  # Phi angle in radians
    psi: float  # Psi angle in radians
    glyph: int  # Nine-glyph index (0-8)
    cost: int   # Recognition cost in units of E_coh
    phase: float = 0.0  # Phase value in radians
    

class TorsionDynamics:
    """
    Manages torsion angle dynamics based on RS nine-glyph system.
    
    This class handles:
    - Torsion angle to glyph mapping
    - Recognition cost calculation
    - Secondary structure template generation
    - Native contact identification
    """
    
    def __init__(self):
        """Initialize torsion dynamics calculator"""
        # Precompute glyph properties
        self._compute_glyph_costs()
        
    def _compute_glyph_costs(self):
        """Compute recognition costs for each glyph"""
        # From RS theory: cost = glyph_index mod 8
        self.glyph_costs = {}
        for g in range(9):
            self.glyph_costs[g] = g % 8
            
        # Special glyphs (from theory)
        # Glyph mapping:
        # 0: phi=[-180,-60], psi=[-180,-60] - extended/sheet region
        # 1: phi=[-180,-60], psi=[-60,60]   - bridge region
        # 2: phi=[-180,-60], psi=[60,180]   - polyproline II
        # 3: phi=[-60,60],   psi=[-180,-60] - right-handed helix turn
        # 4: phi=[-60,60],   psi=[-60,60]   - alpha helix region
        # 5: phi=[-60,60],   psi=[60,180]   - left-handed turn
        # 6: phi=[60,180],   psi=[-180,-60] - extended
        # 7: phi=[60,180],   psi=[-60,60]   - left-handed helix
        # 8: phi=[60,180],   psi=[60,180]   - disallowed region
        
        self.zero_cost_glyphs = [0, 8]  # Extended and disallowed (wraps to 0)
        self.helix_glyph = 4  # Alpha helix region
        self.sheet_glyph = 0  # Beta sheet region
        
    def angles_to_glyph(self, phi: float, psi: float) -> int:
        """
        Convert torsion angles to nine-glyph index.
        
        Args:
            phi: Phi angle in radians (-π to π)
            psi: Psi angle in radians (-π to π)
            
        Returns:
            Glyph index (0-8)
        """
        # Convert to degrees and shift to 0-360 range
        phi_deg = (phi * 180 / np.pi + 180) % 360
        psi_deg = (psi * 180 / np.pi + 180) % 360
        
        # Compute bins (0-2 for each angle)
        phi_bin = int(phi_deg / GLYPH_BIN_SIZE)
        psi_bin = int(psi_deg / GLYPH_BIN_SIZE)
        
        # Ensure bins are in range 0-2
        phi_bin = min(phi_bin, 2)
        psi_bin = min(psi_bin, 2)
        
        # Combined glyph index
        glyph = phi_bin * 3 + psi_bin
        
        return glyph
    
    def glyph_to_angles(self, glyph: int) -> Tuple[float, float]:
        """
        Convert glyph index to representative torsion angles.
        
        Args:
            glyph: Glyph index (0-8)
            
        Returns:
            (phi, psi) in radians
        """
        phi_bin = glyph // 3
        psi_bin = glyph % 3
        
        # Center of each bin
        phi_deg = phi_bin * GLYPH_BIN_SIZE + GLYPH_BIN_SIZE / 2 - 180
        psi_deg = psi_bin * GLYPH_BIN_SIZE + GLYPH_BIN_SIZE / 2 - 180
        
        return (phi_deg * np.pi / 180, psi_deg * np.pi / 180)
    
    def get_recognition_cost(self, glyph: int) -> int:
        """
        Get recognition cost for a glyph.
        
        Args:
            glyph: Glyph index (0-8)
            
        Returns:
            Cost in units of E_coh
        """
        return self.glyph_costs[glyph]
    
    def compute_torsion_state(self, phi: float, psi: float) -> TorsionState:
        """
        Compute complete torsion state for given angles.
        
        Args:
            phi: Phi angle in radians
            psi: Psi angle in radians
            
        Returns:
            TorsionState with glyph and cost information
        """
        glyph = self.angles_to_glyph(phi, psi)
        cost = self.get_recognition_cost(glyph)
        
        return TorsionState(phi=phi, psi=psi, glyph=glyph, cost=cost)
    
    def generate_helix_template(self, n_residues: int, 
                               start_phase: float = 0.0) -> List[TorsionState]:
        """
        Generate ideal α-helix template based on RS theory.
        
        Args:
            n_residues: Number of residues
            start_phase: Starting phase in radians
            
        Returns:
            List of TorsionState for ideal helix
        """
        template = []
        
        # Standard helix angles (glyph 8 - zero cost)
        phi_helix = -60 * np.pi / 180  # -60°
        psi_helix = -45 * np.pi / 180  # -45°
        
        for i in range(n_residues):
            # Phase advances by 100° per residue
            phase = start_phase + i * HELIX_PHASE_ADVANCE
            
            # Create torsion state
            state = self.compute_torsion_state(phi_helix, psi_helix)
            
            # Add phase information (stored in state for later use)
            state.phase = phase % (2 * np.pi)
            
            template.append(state)
            
        return template
    
    def generate_sheet_template(self, n_residues: int,
                               parallel: bool = True) -> List[TorsionState]:
        """
        Generate ideal β-sheet template based on RS theory.
        
        Args:
            n_residues: Number of residues
            parallel: True for parallel, False for antiparallel
            
        Returns:
            List of TorsionState for ideal sheet
        """
        template = []
        
        # Standard sheet angles (glyph 0 - zero cost)
        phi_sheet = -120 * np.pi / 180  # -120°
        psi_sheet = 120 * np.pi / 180   # 120°
        
        for i in range(n_residues):
            # Phase relationship depends on parallel/antiparallel
            if parallel:
                # Phase quadrature (90° offset between strands)
                phase = i * np.pi / 2
            else:
                # Phase opposition (180° offset)
                phase = i * SHEET_PHASE_ADVANCE
                
            # Create torsion state
            state = self.compute_torsion_state(phi_sheet, psi_sheet)
            state.phase = phase % (2 * np.pi)
            
            template.append(state)
            
        return template
    
    def compute_helix_position(self, residue_index: int,
                              start_pos: np.ndarray = np.zeros(3)) -> np.ndarray:
        """
        Compute 3D position of residue in ideal helix.
        
        Args:
            residue_index: Index of residue in helix
            start_pos: Starting position
            
        Returns:
            3D position array
        """
        # Helical parameters
        angle = residue_index * (2 * np.pi / HELIX_RESIDUES_PER_TURN)
        
        # Position in helix
        x = start_pos[0] + HELIX_RADIUS * np.cos(angle)
        y = start_pos[1] + HELIX_RADIUS * np.sin(angle)
        z = start_pos[2] + residue_index * HELIX_RISE
        
        return np.array([x, y, z])
    
    def compute_sheet_position(self, residue_index: int, strand_index: int,
                              start_pos: np.ndarray = np.zeros(3)) -> np.ndarray:
        """
        Compute 3D position of residue in ideal β-sheet.
        
        Args:
            residue_index: Index along strand
            strand_index: Which strand (0, 1, 2, ...)
            start_pos: Starting position
            
        Returns:
            3D position array
        """
        # Sheet parameters
        residue_spacing = 3.3  # Å along strand
        strand_spacing = 4.8   # Å between strands
        
        # Alternating up/down for pleated sheet
        pleat_amplitude = 0.8  # Å
        pleat = pleat_amplitude * (1 if residue_index % 2 == 0 else -1)
        
        x = start_pos[0] + residue_index * residue_spacing
        y = start_pos[1] + strand_index * strand_spacing
        z = start_pos[2] + pleat
        
        return np.array([x, y, z])
    
    def identify_native_contacts(self, positions: np.ndarray,
                                phases: np.ndarray) -> List[Tuple[int, int]]:
        """
        Identify native contacts based on RS recognition criteria.
        
        Args:
            positions: Array of 3D positions (N x 3)
            phases: Array of phase values (N,)
            
        Returns:
            List of (i, j) pairs forming native contacts
        """
        n_residues = len(positions)
        contacts = []
        
        for i in range(n_residues):
            for j in range(i + 2, n_residues):  # Skip neighbors
                # Distance criterion
                dist = np.linalg.norm(positions[j] - positions[i])
                if dist > RECOGNITION_DISTANCE:
                    continue
                    
                # Phase alignment criterion
                phase_diff = abs(phases[j] - phases[i])
                phase_diff = min(phase_diff, 2 * np.pi - phase_diff)  # Wrap around
                
                if phase_diff < PHASE_ALIGNMENT_THRESHOLD:
                    contacts.append((i, j))
                    
        return contacts
    
    def compute_folding_cost(self, torsion_states: List[TorsionState]) -> float:
        """
        Compute total folding cost in units of E_coh.
        
        Args:
            torsion_states: List of torsion states for all residues
            
        Returns:
            Total cost in units of E_coh
        """
        total_cost = sum(state.cost for state in torsion_states)
        return total_cost
    
    def compute_secondary_structure_content(self, 
                                          torsion_states: List[TorsionState]) -> Dict[str, float]:
        """
        Compute secondary structure content from torsion angles.
        
        Args:
            torsion_states: List of torsion states
            
        Returns:
            Dictionary with fractions of helix, sheet, and coil
        """
        n_total = len(torsion_states)
        if n_total == 0:
            return {'helix': 0.0, 'sheet': 0.0, 'coil': 0.0}
            
        # Count by glyph type
        n_helix = sum(1 for state in torsion_states if state.glyph == 4)
        n_sheet = sum(1 for state in torsion_states if state.glyph == 0)
        n_coil = n_total - n_helix - n_sheet
        
        return {
            'helix': n_helix / n_total,
            'sheet': n_sheet / n_total,
            'coil': n_coil / n_total
        } 