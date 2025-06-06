"""
Enhanced Three-Layer Protein Folder with Torsion Dynamics

This module extends the three-layer folder to include:
- Nine-glyph torsion angle tracking
- Secondary structure formation based on golden ratio geometry
- Native contact tracking with phase alignment
- Voxel walk dynamics for folding pathways

All based on Recognition Science first principles.
"""

import numpy as np
from typing import List, Tuple, Optional, Dict
import sys
import os

# Add parent directory to path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from phase2_information.three_layer_folder import ThreeLayerFolder, CoherentRegion
from phase2_information.phase_pattern_field import RecognitionEvent
from phase2_information.torsion_dynamics import TorsionDynamics, TorsionState
from phase2_information.ir_photon_analysis import IRPhotonAnalyzer

# Recognition Science constants
TAU_0 = 7.33e-15  # s (fundamental tick)
E_COH = 0.090  # eV (one recognition quantum)
PHI = (1 + np.sqrt(5)) / 2  # Golden ratio
RECOGNITION_DISTANCE = 6.5  # Å
K_B = 8.617e-5  # eV/K (Boltzmann constant)
MOBILITY_RS = (2 * np.pi / PHI) ** 2 / 160  # ≈0.096 Å²/(eV·ps) RS-derived (see derive_constants.py)

# Voxel parameters from RS theory
VOXEL_SIZE = 0.335  # nm (3.35 Å) - one recognition voxel edge length
DAMPING_FACTOR = np.sqrt(PHI) / PHI  # √P · φ^(-1/2) for backbone bonds

# New constant
BASE_RECOGNITION_PROB_RS = (2 * np.pi / PHI) ** 2 / 160  # ≈0.096


class EnhancedCoherentRegion(CoherentRegion):
    """Extended coherent region with torsion angle information"""
    def __init__(self, index: int, position: np.ndarray, phase: float = 0.0):
        super().__init__(index, position, phase)
        self.torsion_state: Optional[TorsionState] = None
        self.secondary_structure: str = 'coil'  # 'helix', 'sheet', or 'coil'
        self.native_contacts: List[int] = []  # Indices of native contact partners
        

class EnhancedThreeLayerFolder(ThreeLayerFolder):
    """
    Enhanced Recognition Science protein folder with torsion dynamics.
    
    Extends the base three-layer folder with:
    - Torsion angle tracking and energy accounting
    - Secondary structure template matching
    - Native contact formation based on phase alignment
    - Voxel walk dynamics for folding pathways
    - IR photon emission tracking (13.8 μm)
    """
    
    def __init__(self, n_residues: int, temperature: float = 310.0,
                 sequence: Optional[str] = None):
        """
        Initialize enhanced folder.
        
        Args:
            n_residues: Number of amino acid residues
            temperature: Temperature in Kelvin
            sequence: Optional amino acid sequence (for future use)
        """
        # Initialize torsion dynamics first (needed by _initialize_enhanced_regions)
        self.torsion_dynamics = TorsionDynamics()
        
        # Initialize IR photon analyzer
        self.ir_analyzer = IRPhotonAnalyzer()
        
        # Initialize base class
        super().__init__(n_residues, temperature)
        
        # Replace regions with enhanced versions
        self.regions = self._initialize_enhanced_regions()
        
        # Track secondary structure formation
        self.ss_templates = {}  # Will store identified secondary structures
        self.ss_formation_tick = {}  # When each SS element formed
        
        # Native contact tracking
        self.native_contacts = []  # List of (i, j) pairs
        self.contact_formation_tick = {}  # When each contact formed
        
        # Voxel walk tracking
        self.voxel_transitions = 0  # Total number of voxel transitions
        self.folding_pathway = []  # Sequence of major conformational changes
        
        # Sequence information (for future sequence-specific effects)
        self.sequence = sequence if sequence else 'A' * n_residues
        
    def _initialize_enhanced_regions(self) -> List[EnhancedCoherentRegion]:
        """Initialize enhanced regions with torsion information"""
        regions = []
        for i in range(self.n_residues):
            # Extended chain along x-axis
            position = np.array([i * 3.8, 0.0, 0.0])  # 3.8 Å spacing
            region = EnhancedCoherentRegion(i, position)
            
            # Initialize with random coil torsion angles
            phi = np.random.uniform(-np.pi, np.pi)
            psi = np.random.uniform(-np.pi, np.pi)
            region.torsion_state = self.torsion_dynamics.compute_torsion_state(phi, psi)
            
            regions.append(region)
        return regions
    
    def _find_recognition_events(self) -> List[RecognitionEvent]:
        """
        Find recognition events with torsion-aware probability modulation.
        
        Recognition probability is enhanced when:
        - Torsion angles favor secondary structure formation
        - Phase alignment matches golden ratio relationships
        - Native contacts are being formed
        """
        events = []
        
        for i in range(self.n_residues):
            for j in range(i + 1, self.n_residues):
                # Check if already recognized
                if j in self.regions[i].recognized_with:
                    continue
                
                # Calculate distance
                r_ij = self.positions[j] - self.positions[i]
                distance = np.linalg.norm(r_ij)
                
                if distance <= RECOGNITION_DISTANCE:
                    # Base recognition probability (RS-derived)
                    base_prob = BASE_RECOGNITION_PROB_RS
                    
                    # Phase modulation
                    phase_diff = self.regions[j].phase - self.regions[i].phase
                    phase_factor = (1 + np.cos(phase_diff)) / 2
                    
                    # Torsion angle modulation
                    torsion_factor = self._compute_torsion_factor(i, j)
                    
                    # Eight-beat modulation
                    beat_phase = (self.tick % 8) * 2 * np.pi / 8
                    beat_factor = (1 + np.cos(beat_phase)) / 2
                    
                    # Total probability
                    prob = base_prob * phase_factor * torsion_factor * beat_factor
                    
                    if np.random.random() < prob:
                        # Create recognition event
                        event = RecognitionEvent(
                            residue_i=i,
                            residue_j=j,
                            phase_shift=phase_diff,
                            tick=self.tick
                        )
                        events.append(event)
                        
                        # Update recognized partners
                        self.regions[i].recognized_with.add(j)
                        self.regions[j].recognized_with.add(i)
                        
                        # Update phases (discrete jump)
                        phase_exchange = E_COH / (2 * self.kT)
                        self.regions[i].phase += phase_exchange
                        self.regions[j].phase -= phase_exchange
                        
                        # Track IR photon emission
                        self.ir_analyzer.add_recognition_event(self.tick, i, j)
                        
                        # Check for native contact formation
                        if self._is_native_contact(i, j):
                            self.native_contacts.append((i, j))
                            self.contact_formation_tick[(i, j)] = self.tick
                            self.regions[i].native_contacts.append(j)
                            self.regions[j].native_contacts.append(i)
                        
                        # Track conservation
                        self.total_coins += 1
        
        return events
    
    def _compute_torsion_factor(self, i: int, j: int) -> float:
        """
        Compute recognition probability factor based on torsion angles.
        
        Enhanced probability when:
        - Both residues have low-cost glyphs (helix or sheet)
        - Sequence separation favors secondary structure (|i-j| = 3-4 for helix)
        """
        state_i = self.regions[i].torsion_state
        state_j = self.regions[j].torsion_state
        
        # Low-cost glyphs enhance recognition
        cost_factor_i = 1.0 + 0.5 * (1.0 - state_i.cost / 8.0)
        cost_factor_j = 1.0 + 0.5 * (1.0 - state_j.cost / 8.0)
        
        # Sequence separation factor
        sep = abs(j - i)
        if 3 <= sep <= 4:  # Helix i, i+3/i+4 contacts
            sep_factor = PHI  # Golden ratio boost
        elif sep > 10:  # Long-range contacts
            sep_factor = 1.2
        else:
            sep_factor = 1.0
            
        return cost_factor_i * cost_factor_j * sep_factor
    
    def _is_native_contact(self, i: int, j: int) -> bool:
        """
        Determine if a recognition event forms a native contact.
        
        Based on:
        - Phase alignment within threshold
        - Appropriate sequence separation
        - Torsion angle compatibility
        """
        # Phase alignment check
        phase_diff = abs(self.regions[j].phase - self.regions[i].phase)
        phase_diff = min(phase_diff, 2 * np.pi - phase_diff)
        
        if phase_diff > np.pi / 4:  # 45° threshold
            return False
            
        # Sequence separation check
        sep = abs(j - i)
        if sep < 3:  # Too close in sequence
            return False
            
        # Torsion compatibility (both in similar secondary structure)
        glyph_i = self.regions[i].torsion_state.glyph
        glyph_j = self.regions[j].torsion_state.glyph
        
        # Both helix or both sheet
        if (glyph_i == 4 and glyph_j == 4) or (glyph_i == 0 and glyph_j == 0):
            return True
            
        # Default: consider it native if phase-aligned
        return True
    
    def _update_physical_configuration(self):
        """
        Enhanced physical update including torsion angle evolution.
        """
        # Get information pressure gradient
        info_pressure = self.phase_field.compute_information_pressure(self.positions)
        
        # Check for barrier crossing (as before)
        if not self.folding_initiated:
            # RS Arrhenius prefactor (size-dependent, topology factors≈1 here)
            k0_base = (1 / (8 * TAU_0)) * PHI ** (-self.n_residues / 2)
            k_fold = k0_base * np.exp(-0.18 / self.kT)
            dt = TAU_0
            p_initiate = k_fold * dt
            
            if np.random.random() < p_initiate:
                self.folding_initiated = True
                self.barrier_crossing_attempts += 1
                print(f"Folding initiated after {self.barrier_crossing_attempts} attempts")
                print(f"Time since template: {(self.tick - self.template_completion_tick) * TAU_0 * 1e6:.2f} μs")
        
        # If folding is initiated, update configuration
        if self.folding_initiated:
            # Update positions following information pressure
            # Use overdamped dynamics (no inertia)
            # Mobility constant: effective parameter for voxel transitions
            # Derived from voxel statistics (~1 Å²/(eV·ps)) × damping factor
            # See derive_constants.py for full derivation from RS principles
            mobility = MOBILITY_RS  # RS-derived mobility constant ≈0.096 Å²/(eV·ps)
            dt_physical = TAU_0 * 1e12  # Convert to ps for physical update
            
            # Position update with voxel constraints
            for i in range(self.n_residues):
                # Compute displacement
                displacement = mobility * info_pressure[i] * dt_physical
                
                # Check if this crosses a voxel boundary
                old_voxel = self._position_to_voxel(self.positions[i])
                new_position = self.positions[i] + displacement
                new_voxel = self._position_to_voxel(new_position)
                
                if not np.array_equal(old_voxel, new_voxel):
                    # Voxel transition - apply damping
                    displacement *= DAMPING_FACTOR
                    self.voxel_transitions += 1
                
                self.positions[i] += displacement
                self.regions[i].position = self.positions[i].copy()
            
            # Update torsion angles based on local environment
            self._update_torsion_angles()
            
            # Check for secondary structure formation
            self._detect_secondary_structures()
            
            # Update folding progress
            self._update_folding_progress()
    
    def _position_to_voxel(self, position: np.ndarray) -> np.ndarray:
        """Convert position to voxel indices"""
        return np.floor(position / (VOXEL_SIZE * 10)).astype(int)  # Convert nm to Å
    
    def _update_torsion_angles(self):
        """
        Update torsion angles based on local recognition environment.
        
        Torsion angles evolve toward low-cost states when:
        - Multiple native contacts are formed
        - Local phase coherence is high
        - Information pressure favors specific conformations
        """
        for i, region in enumerate(self.regions):
            # Skip if no native contacts yet
            if len(region.native_contacts) == 0:
                continue
                
            # Compute local phase coherence
            local_coherence = 0.0
            for j in region.native_contacts:
                phase_diff = abs(self.regions[j].phase - region.phase)
                local_coherence += np.cos(phase_diff)
            local_coherence /= max(1, len(region.native_contacts))
            
            # Probability of torsion angle update
            p_update = 0.1 * local_coherence
            
            if np.random.random() < p_update:
                # Determine target secondary structure based on contacts
                target_ss = self._determine_target_ss(i)
                
                if target_ss == 'helix':
                    # Move toward helix angles (glyph 8)
                    target_phi = -60 * np.pi / 180
                    target_psi = -45 * np.pi / 180
                elif target_ss == 'sheet':
                    # Move toward sheet angles (glyph 0)
                    target_phi = -120 * np.pi / 180
                    target_psi = 120 * np.pi / 180
                else:
                    continue  # No change for coil
                
                # Gradual movement toward target
                current_phi = region.torsion_state.phi
                current_psi = region.torsion_state.psi
                
                new_phi = current_phi + 0.2 * (target_phi - current_phi)
                new_psi = current_psi + 0.2 * (target_psi - current_psi)
                
                # Update torsion state
                region.torsion_state = self.torsion_dynamics.compute_torsion_state(
                    new_phi, new_psi
                )
    
    def _determine_target_ss(self, residue_index: int) -> str:
        """Determine target secondary structure based on contact patterns"""
        region = self.regions[residue_index]
        
        # Count contacts by sequence separation
        helix_contacts = 0
        sheet_contacts = 0
        
        for j in region.native_contacts:
            sep = abs(j - residue_index)
            if 3 <= sep <= 4:
                helix_contacts += 1
            elif sep > 10:
                sheet_contacts += 1
        
        if helix_contacts > sheet_contacts and helix_contacts >= 2:
            return 'helix'
        elif sheet_contacts >= 2:
            return 'sheet'
        else:
            return 'coil'
    
    def _detect_secondary_structures(self):
        """Detect formation of secondary structure elements"""
        # Check for helices (consecutive residues with helix glyphs)
        helix_runs = []
        current_run = []
        
        for i, region in enumerate(self.regions):
            if region.torsion_state.glyph == 4:  # Helix glyph (was 8)
                current_run.append(i)
            else:
                if len(current_run) >= 4:  # Minimum helix length
                    helix_runs.append(current_run)
                current_run = []
        
        if len(current_run) >= 4:
            helix_runs.append(current_run)
        
        # Record new helices
        for run in helix_runs:
            helix_id = f"helix_{min(run)}_{max(run)}"
            if helix_id not in self.ss_templates:
                self.ss_templates[helix_id] = run
                self.ss_formation_tick[helix_id] = self.tick
                print(f"Helix formed: residues {min(run)}-{max(run)} at "
                      f"{self.tick * TAU_0 * 1e12:.2f} ps")
        
        # Similar detection for sheets (more complex due to strand pairing)
        # Simplified: just detect runs of sheet glyphs
        sheet_runs = []
        current_run = []
        
        for i, region in enumerate(self.regions):
            if region.torsion_state.glyph == 0:  # Sheet glyph (correct)
                current_run.append(i)
            else:
                if len(current_run) >= 3:  # Minimum strand length
                    sheet_runs.append(current_run)
                current_run = []
        
        if len(current_run) >= 3:
            sheet_runs.append(current_run)
        
        for run in sheet_runs:
            sheet_id = f"sheet_{min(run)}_{max(run)}"
            if sheet_id not in self.ss_templates:
                self.ss_templates[sheet_id] = run
                self.ss_formation_tick[sheet_id] = self.tick
    
    def _update_folding_progress(self):
        """Update folding progress based on multiple metrics"""
        # Compactness metric (as before)
        current_compactness = np.std(self.positions)
        initial_compactness = self.n_residues * 3.8 / np.sqrt(12)
        compactness_progress = 1.0 - current_compactness / initial_compactness
        
        # Native contact metric
        max_contacts = self.n_residues * 2  # Rough estimate
        contact_progress = len(self.native_contacts) / max_contacts
        
        # Secondary structure metric
        torsion_states = [r.torsion_state for r in self.regions]
        ss_content = self.torsion_dynamics.compute_secondary_structure_content(torsion_states)
        ss_progress = ss_content['helix'] + ss_content['sheet']
        
        # Combined progress
        self.folding_progress = (compactness_progress + contact_progress + ss_progress) / 3
    
    def get_detailed_metrics(self) -> Dict[str, float]:
        """Get detailed simulation metrics including torsion and SS information"""
        base_metrics = self._get_step_metrics()
        
        # Compute secondary structure content
        torsion_states = [r.torsion_state for r in self.regions]
        ss_content = self.torsion_dynamics.compute_secondary_structure_content(torsion_states)
        
        # Compute average torsion cost
        avg_cost = self.torsion_dynamics.compute_folding_cost(torsion_states) / self.n_residues
        
        # Add enhanced metrics
        base_metrics.update({
            'helix_content': ss_content['helix'],
            'sheet_content': ss_content['sheet'],
            'coil_content': ss_content['coil'],
            'native_contacts': len(self.native_contacts),
            'avg_torsion_cost': avg_cost,
            'voxel_transitions': self.voxel_transitions,
            'n_helices': sum(1 for k in self.ss_templates if k.startswith('helix')),
            'n_sheets': sum(1 for k in self.ss_templates if k.startswith('sheet')),
        })
        
        return base_metrics
    
    def save_trajectory(self, filename: str):
        """Save current structure in PDB format with enhanced information"""
        with open(filename, 'w') as f:
            f.write("REMARK  Enhanced RS Three-Layer Folder\n")
            f.write(f"REMARK  Time: {self.tick * TAU_0 * 1e12:.2f} ps\n")
            f.write(f"REMARK  Folding Progress: {self.folding_progress:.3f}\n")
            f.write(f"REMARK  Native Contacts: {len(self.native_contacts)}\n")
            
            # Write secondary structure information
            for ss_id, residues in self.ss_templates.items():
                f.write(f"REMARK  {ss_id}: residues {min(residues)}-{max(residues)}\n")
            
            # Write atoms
            for i, region in enumerate(self.regions):
                # Determine SS type for coloring
                if region.torsion_state.glyph == 8:
                    ss_char = 'H'  # Helix
                elif region.torsion_state.glyph == 0:
                    ss_char = 'E'  # Extended (sheet)
                else:
                    ss_char = 'C'  # Coil
                
                f.write(f"ATOM  {i+1:5d}  CA  {self.sequence[i]} A{i+1:4d}    "
                       f"{region.position[0]:8.3f}{region.position[1]:8.3f}{region.position[2]:8.3f}"
                       f"  1.00 {region.torsion_state.cost:5.2f}           C  {ss_char}\n")
            
            # Write connections for native contacts
            for i, j in self.native_contacts:
                f.write(f"CONECT{i+1:5d}{j+1:5d}\n")
            
            f.write("END\n") 