"""
Pattern Analyzer - Extract protein-specific parameters from formed templates

This module analyzes the PhasePatternField after template formation to extract
topological features that determine folding kinetics. This is the RS-pure approach:
let the emerged pattern tell us the parameters rather than guessing from sequence.
"""

import numpy as np
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass
from collections import defaultdict
import networkx as nx

# RS constants
PHI = 1.618033988749895
E_COH = 0.090  # eV

@dataclass
class TemplateAnalysis:
    """Results from analyzing a formed template."""
    barrier_coins: int
    barrier_ev: float
    p_ledger: float
    p_geom: float
    n_voxels: int
    n_components: int
    n_loops: int
    unique_rungs: Set[int]
    contact_order: float
    helix_fraction: float
    sheet_fraction: float
    path_entropy: float = 1.0
    mobility_anisotropy: float = 1.0
    
class PatternAnalyzer:
    """Analyzes PhasePatternField to extract topological parameters."""
    
    def __init__(self, voxel_size: float = 0.335e-9):
        """
        Initialize analyzer.
        
        Args:
            voxel_size: Size of recognition voxel in meters (default 3.35 Å)
        """
        self.voxel_size = voxel_size
        
    def analyze_template(self, phase_field, positions: np.ndarray, 
                        torsion_states: np.ndarray) -> TemplateAnalysis:
        """
        Extract parameters from a formed template.
        
        Args:
            phase_field: PhasePatternField after template formation
            positions: Current atom positions [n_atoms, 3]
            torsion_states: Current torsion glyph indices [n_residues]
            
        Returns:
            TemplateAnalysis with all extracted parameters
        """
        # Store phase_field for use in methods
        self.phase_field = phase_field
        
        # Build voxel graph from occupied cells
        voxel_graph, voxel_map = self._build_voxel_graph(phase_field, positions)
        
        # Count topological features
        n_components = nx.number_connected_components(voxel_graph)
        n_loops = self._count_independent_loops(voxel_graph)
        
        # Identify unique torsion rungs
        unique_rungs = self._identify_unique_rungs(torsion_states)
        
        # Calculate contact order
        contact_order = self._calculate_contact_order(voxel_graph, voxel_map)
        
        # Detect secondary structure content
        helix_frac, sheet_frac = self._analyze_secondary_structure(torsion_states)
        
        # RS topological counts
        helix_axes = self._count_helix_axes(torsion_states)
        beta_registries = self._count_beta_registries(voxel_graph, voxel_map, torsion_states, positions)
        
        # Debug output
        print(f"DEBUG: Helix axes (L) = {helix_axes}, Beta registries (R) = {beta_registries}")
        print(f"DEBUG: Sheet fraction = {sheet_frac:.2f}")
        print(f"DEBUG: Torsion states = {torsion_states}")
        sheet_residues = [i for i, t in enumerate(torsion_states) if t in [0, 8]]
        print(f"DEBUG: Sheet residues at positions: {sheet_residues}")
        
        # NEW: Path entropy and mobility tensor (stubbed for now)
        path_entropy = self._estimate_path_entropy(voxel_graph)
        mobility_anisotropy = self._estimate_mobility_anisotropy(voxel_graph)
        
        # Calculate parameters from topology
        barrier_coins = self._calculate_barrier(unique_rungs, n_components)
        p_ledger = self._calculate_p_ledger(n_components)
        p_geom = self._calculate_p_geom(n_loops, contact_order, helix_axes, beta_registries)
        
        return TemplateAnalysis(
            barrier_coins=barrier_coins,
            barrier_ev=barrier_coins * E_COH,
            p_ledger=p_ledger,
            p_geom=p_geom,
            n_voxels=voxel_graph.number_of_nodes(),
            n_components=n_components,
            n_loops=n_loops,
            unique_rungs=unique_rungs,
            contact_order=contact_order,
            helix_fraction=helix_frac,
            sheet_fraction=sheet_frac,
            path_entropy=path_entropy,
            mobility_anisotropy=mobility_anisotropy
        )
    
    def _build_voxel_graph(self, phase_field, positions: np.ndarray) -> Tuple[nx.Graph, Dict]:
        """
        Build graph of voxel connectivity from phase field.
        
        Returns:
            (graph, voxel_map) where voxel_map[voxel_idx] = residue_indices
        """
        # Get voxel indices for each atom
        voxel_indices = phase_field.get_voxel_indices(positions)
        
        # Group atoms by voxel
        voxel_map = defaultdict(list)
        for i, vox_idx in enumerate(voxel_indices):
            if vox_idx is not None:  # Skip atoms outside grid
                voxel_map[vox_idx].append(i)
        
        # Build graph where nodes are occupied voxels
        G = nx.Graph()
        occupied_voxels = list(voxel_map.keys())
        G.add_nodes_from(occupied_voxels)
        
        # 1) Add edges between spatially adjacent voxels (26-connectivity)
        for i, vox1 in enumerate(occupied_voxels):
            for vox2 in occupied_voxels[i+1:]:
                if self._are_adjacent_voxels(vox1, vox2, phase_field.grid_size):
                    G.add_edge(vox1, vox2)
        
        # 2) Add edges implied by recognition (phase-reach rule)
        #    Two voxels are connected if any residue pair spanning those voxels
        #    has already exchanged recognition quanta (ledger reach ≥ φ⁻¹).
        # Build residue→voxel map (use first encountered atom)
        residue_to_vox = {}
        for vox_idx, atom_list in voxel_map.items():
            for atom_idx in atom_list:
                res_idx = atom_idx // 4  # four atoms per residue (approx.)
                if res_idx not in residue_to_vox:
                    residue_to_vox[res_idx] = vox_idx
        
        # Iterate over recognised residue pairs
        if hasattr(phase_field, 'recognition_matrix'):
            recog = phase_field.recognition_matrix
            n_residues = recog.shape[0]
            for i in range(n_residues):
                for j in range(i+1, n_residues):
                    if recog[i, j]:
                        vi = residue_to_vox.get(i)
                        vj = residue_to_vox.get(j)
                        if vi is not None and vj is not None and vi != vj:
                            G.add_edge(vi, vj)
        
        return G, dict(voxel_map)
    
    def _are_adjacent_voxels(self, vox1: Tuple, vox2: Tuple, grid_size: np.ndarray) -> bool:
        """Check if two voxel indices are adjacent (26-connectivity)."""
        diff = np.abs(np.array(vox1) - np.array(vox2))
        return np.all(diff <= 1) and np.sum(diff) > 0
    
    def _count_independent_loops(self, graph: nx.Graph) -> int:
        """Count independent loops using cycle basis."""
        if graph.number_of_edges() == 0:
            return 0
        
        # For a connected graph: loops = edges - nodes + 1
        # For multiple components, sum over each
        n_loops = 0
        for component in nx.connected_components(graph):
            subgraph = graph.subgraph(component)
            n_loops += subgraph.number_of_edges() - subgraph.number_of_nodes() + 1
        
        return n_loops
    
    def _identify_unique_rungs(self, torsion_states: np.ndarray) -> Set[int]:
        """Identify unique torsion rungs present."""
        # Torsion states are glyph indices 0-8
        # Map to rungs based on recognition cost
        rungs = set()
        for state in torsion_states:
            if 0 <= state <= 8:
                # Recognition cost = glyph_index mod 8
                rung = state % 8
                rungs.add(rung)
        return rungs
    
    def _calculate_contact_order(self, graph: nx.Graph, voxel_map: Dict) -> float:
        """
        Calculate relative contact order from voxel graph.
        
        Contact order = average sequence separation of contacts / total length
        """
        if graph.number_of_edges() == 0:
            return 0.0
        
        total_separation = 0
        n_contacts = 0
        
        # For each edge in voxel graph
        for vox1, vox2 in graph.edges():
            # Get residue indices in each voxel
            res1 = voxel_map.get(vox1, [])
            res2 = voxel_map.get(vox2, [])
            
            # Calculate minimum sequence separation
            if res1 and res2:
                # Assuming residue index ~ atom index / 4 (rough)
                min_sep = float('inf')
                for r1 in res1:
                    for r2 in res2:
                        sep = abs(r1//4 - r2//4)
                        if sep > 2:  # Skip local contacts
                            min_sep = min(min_sep, sep)
                
                if min_sep < float('inf'):
                    total_separation += min_sep
                    n_contacts += 1
        
        if n_contacts == 0:
            return 0.0
        
        # Estimate total length from number of voxels
        n_residues = len(voxel_map)
        return total_separation / (n_contacts * n_residues)
    
    def _analyze_secondary_structure(self, torsion_states: np.ndarray) -> Tuple[float, float]:
        """
        Detect secondary structure from torsion patterns.
        
        Returns:
            (helix_fraction, sheet_fraction)
        """
        n_residues = len(torsion_states)
        if n_residues < 3:
            return 0.0, 0.0
        
        helix_count = 0
        sheet_count = 0
        
        # Helix: glyphs 2,3 (φ ≈ -60°, ψ ≈ -45°)
        # Sheet: glyphs 0,8 (φ ≈ -120°, ψ ≈ 120°)
        
        for i in range(n_residues):
            if torsion_states[i] in [2, 3]:
                helix_count += 1
            elif torsion_states[i] in [0, 8]:
                sheet_count += 1
        
        return helix_count / n_residues, sheet_count / n_residues
    
    def _calculate_barrier(self, unique_rungs: Set[int], n_components: int) -> int:
        """
        Return the universal Recognition-Science folding barrier:
        ΔC = 2 × E_COH (two coins).  RS axioms forbid protein-specific
        barrier inflation; size/topology differences only affect the
        Arrhenius pre-factor, not ΔC itself.
        """
        return 2
    
    def _calculate_p_ledger(self, n_components: int) -> float:
        """
        Calculate ledger availability factor.
        
        Multiple components require coordinated ledger access.
        """
        # Each component reduces availability by φ^(-1/2)
        return PHI ** (-n_components / 2)
    
    def _calculate_p_geom(self, n_loops: int, contact_order: float,
                         n_helix_axes: int, n_beta_regs: int) -> float:
        """
        Calculate geometric compatibility factor.
        
        Considers:
        - Loop complexity
        - Long-range contacts
        - Secondary structure bonuses/penalties
        """
        # RS-derived formula (Protein-Full-2 §15):
        #   P_geom = φ^{-n_loops/2} · φ^{-(√φ)·CO} · φ^{helix_frac/2} · φ^{-sheet_frac}

        p_geom = PHI ** (-n_loops / 2)

        # Contact-order penalty scales with √φ
        if contact_order > 0:
            p_geom *= PHI ** (-np.sqrt(PHI) * contact_order)

        # Helix-axis bonus (L)
        if n_helix_axes > 0:
            p_geom *= PHI ** (n_helix_axes / 2)

        # Beta-strand registry penalty (R)
        if n_beta_regs > 0:
            p_geom *= PHI ** (-n_beta_regs)

        # RS logic: geometric factor cannot exceed 1 (fully favorable); no
        # empirical 0.1 cap.  Clip only to the theoretical maximum = 1.0.
        if p_geom > 1.0:
            p_geom = 1.0
        return max(p_geom, 1e-8)
    
    # ------------------------------------------------------------------
    #  NEW STUB ESTIMATES (to be replaced with true RS derivations)
    # ------------------------------------------------------------------

    def _estimate_path_entropy(self, graph: nx.Graph) -> float:
        """RS-derived path entropy.

        Definition (see Deeper Understanding.txt §4.3):
            S_path = (1 / N_pairs) * Σ_{u<v} ln[ N_sp(u,v) ]
        where N_sp(u,v) is the number of *distinct* shortest voxel-graph
        paths between voxels u and v.

        We normalise by dividing by ln(φ) so that an ideal maximally
        coherent template (all voxel pairs linked by a *unique* shortest
        path) yields S_path = 0, while highly redundant connectivity
        approaches 1.
        """

        n_nodes = graph.number_of_nodes()
        if n_nodes < 3:
            return 0.0  # Trivial templates

        # Pre-compute all-pairs shortest-path lengths
        sp_lengths = dict(nx.all_pairs_shortest_path_length(graph))

        # Count shortest-path degeneracy for each pair using Yen's algo
        # For small graphs this brute-force approach is acceptable.
        total_log_paths = 0.0
        n_pairs = 0

        for i_idx, u in enumerate(graph.nodes()):
            for v in list(graph.nodes())[i_idx+1:]:
                d = sp_lengths[u].get(v)
                if d is None:
                    continue  # different components
                # Count number of shortest paths of length d
                paths = list(nx.all_shortest_paths(graph, u, v))
                n_paths = max(1, len(paths))
                total_log_paths += np.log(n_paths)
                n_pairs += 1

        if n_pairs == 0:
            return 0.0

        # Normalised entropy (divide by ln(φ) ~ golden-ratio natural log)
        mean_log = total_log_paths / n_pairs
        S_norm = mean_log / np.log(PHI)

        # Clamp to [0,1]
        return float(np.clip(S_norm, 0.0, 1.0))

    def _estimate_mobility_anisotropy(self, graph: nx.Graph) -> float:
        """RS-derived mobility anisotropy.

        We approximate mobility by the eigenspectrum of the voxel cloud
        inertia tensor (covariance of voxel coordinates).  An isotropic
        template has equal eigenvalues → anisotropy 0.

        Anisotropy = (λ_max / λ_min) − 1 (≥0)
        """

        if graph.number_of_nodes() < 3:
            return 0.0

        # Convert voxel indices to coordinates (Å) centred at origin
        coords = np.array(graph.nodes()) * self.voxel_size * 1e10  # Å
        coords = coords - coords.mean(axis=0, keepdims=True)

        # Covariance / inertia tensor
        cov = np.cov(coords.T)
        eigvals = np.linalg.eigvalsh(cov)
        eigvals = np.sort(eigvals)
        if eigvals[0] <= 0:
            return 0.0

        anisotropy = (eigvals[-1] / eigvals[0]) - 1.0
        return float(anisotropy)

    # ------------------------------------------------------------------
    #  Helix-axis (L) and β-registry (R) counters
    # ------------------------------------------------------------------

    def _count_helix_axes(self, torsion_states: np.ndarray) -> int:
        """Count contiguous helix segments ≥6 residues."""
        L = 0
        run_len = 0
        helix_segments = []
        current_start = None
        
        for state in torsion_states:
            if state in [2, 3]:
                if run_len == 0:
                    current_start = len(helix_segments)
                run_len += 1
            else:
                if run_len >= 6:
                    L += 1
                    helix_segments.append(run_len)
                elif run_len > 0:
                    helix_segments.append(f"({run_len})")  # Too short
                run_len = 0
                current_start = None
                
        if run_len >= 6:
            L += 1
            helix_segments.append(run_len)
        elif run_len > 0:
            helix_segments.append(f"({run_len})")
            
        helix_residues = [i for i, t in enumerate(torsion_states) if t in [2, 3]]
        print(f"DEBUG helix: Found {L} helix axes (≥6), segments: {helix_segments}")
        print(f"DEBUG helix: Helix residues at: {helix_residues}")
        
        return L

    def _count_beta_registries(self, graph: nx.Graph, voxel_map: Dict,
                               torsion_states: np.ndarray, positions: np.ndarray) -> int:
        """Count independent β-strand registry matches (R) using recognition matrix.
        
        RS-based algorithm:
        1. Identify contiguous sheet segments ≥3 residues (torsion 0 or 8)
        2. For pairs of segments with sequence separation >2:
           - Check if residues have recognized each other (recognition_matrix[i,j]=True)
           - A ladder exists if ≥2 consecutive residue pairs have recognized
        3. Each ladder pattern counts as one registry (R)
        
        This directly counts recognition events rather than geometric distances,
        staying true to RS first principles.
        """
        # Get recognition matrix from phase field
        if not hasattr(self.phase_field, 'recognition_matrix'):
            return 0
        
        recog_matrix = self.phase_field.recognition_matrix
        n_residues = len(torsion_states)
        
        # Debug: Check recognition matrix
        n_recognitions = np.sum(recog_matrix) // 2  # Symmetric matrix
        print(f"DEBUG β-registry: Recognition matrix has {n_recognitions} recognition events")
        print(f"DEBUG β-registry: Matrix shape = {recog_matrix.shape}")
        
        # Show recognition pairs
        recog_pairs = []
        for i in range(n_residues):
            for j in range(i+1, n_residues):
                if recog_matrix[i, j]:
                    recog_pairs.append((i, j))
        print(f"DEBUG β-registry: Recognition pairs: {recog_pairs[:10]}...")  # First 10
        
        # Step 1: Find sheet segments
        segments = []  # list of (start_idx, end_idx)
        i = 0
        while i < len(torsion_states):
            if torsion_states[i] in [0, 8]:
                start = i
                while i < len(torsion_states) and torsion_states[i] in [0, 8]:
                    i += 1
                end = i - 1
                if end - start + 1 >= 3:  # Minimum 3 residues for a sheet segment
                    segments.append((start, end))
            else:
                i += 1
        
        print(f"DEBUG β-registry: Found {len(segments)} sheet segments: {segments}")
        
        # Count registries from segments
        n_registries = 0
        
        if len(segments) >= 2:
            # Check segment pairs for ladder patterns
            for idx1 in range(len(segments)):
                s1_start, s1_end = segments[idx1]
                for idx2 in range(idx1 + 1, len(segments)):
                    s2_start, s2_end = segments[idx2]
                    
                    # Check sequence separation
                    if s2_start - s1_end <= 2:
                        continue
                    
                    # Look for ladder pattern
                    ladder_found = False
                    
                    # Try all possible alignments
                    for offset in range(-(s2_end - s2_start), (s1_end - s1_start) + 1):
                        consecutive = 0
                        
                        for i in range(max(0, -offset), 
                                     min(s1_end - s1_start + 1, s2_end - s2_start + 1 - offset)):
                            res1 = s1_start + i
                            res2 = s2_start + i + offset
                            
                            if 0 <= res2 < n_residues and recog_matrix[res1, res2]:
                                consecutive += 1
                                if consecutive >= 2:
                                    ladder_found = True
                                    break
                            else:
                                consecutive = 0
                        
                        if ladder_found:
                            break
                    
                    if ladder_found:
                        n_registries += 1
        
        # Alternative: Look for recognition patterns even without clear segments
        # This catches β-hairpins and other patterns missed by torsion analysis
        ladders_found = 0
        
        # Look for anti-parallel ladder patterns
        for i in range(n_residues - 5):  # Need space for ladder
            for j in range(i + 4, n_residues):  # Skip nearby residues
                # Check for at least 2 consecutive recognitions
                if j >= 1 and recog_matrix[i, j] and recog_matrix[i+1, j-1]:
                    ladders_found += 1
                    print(f"DEBUG: Found anti-parallel ladder at {i}-{j}")
                
                # Check for parallel ladder
                elif j < n_residues - 1 and recog_matrix[i, j] and recog_matrix[i+1, j+1]:
                    ladders_found += 1
                    print(f"DEBUG: Found parallel ladder at {i}-{j}")
        
        total_registries = n_registries + ladders_found
        print(f"DEBUG: Total registries = {total_registries} (segments: {n_registries}, patterns: {ladders_found})")
        
        # If no registries found but we have appreciable sheet content (> φ⁻² ~ 0.382),
        # estimate a minimal registry count.  This RS-derived threshold avoids
        # the previous empirical 0.1 value.
        if total_registries == 0 and self._analyze_secondary_structure(torsion_states)[1] > PHI**-2:
            # Count potential β-hairpin sites (common in small proteins)
            # Look for patterns like: helix-turn-sheet or sheet-turn-sheet
            estimated_registries = 0
            
            # Simple heuristic: if we have scattered sheet residues,
            # assume they'll form at least one registry during folding
            sheet_positions = [i for i, t in enumerate(torsion_states) if t in [0, 8]]
            if len(sheet_positions) >= 2:
                # Check if they're separated enough to form a hairpin
                for i in range(len(sheet_positions)-1):
                    if sheet_positions[i+1] - sheet_positions[i] >= 3:
                        estimated_registries = 1
                        break
            
            if estimated_registries > 0:
                print(f"DEBUG: Estimated {estimated_registries} registries from sequence (template too early)")
                return estimated_registries
        
        return total_registries 