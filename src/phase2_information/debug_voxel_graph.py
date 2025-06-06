#!/usr/bin/env python3
"""Debug voxel graph construction in pattern analyzer."""

import numpy as np
from phase_pattern_field import PhasePatternField
from pattern_analyzer import PatternAnalyzer

# Create a simple test case
n_residues = 5
positions = np.zeros((n_residues * 4, 3))  # 4 atoms per residue

# Place residues in a line, 3.5 Å apart
for i in range(n_residues):
    for j in range(4):  # 4 atoms per residue
        atom_idx = i * 4 + j
        positions[atom_idx] = [i * 3.5, j * 1.0, 0.0]  # Å

# Create phase field
phase_field = PhasePatternField(grid_size=(10, 10, 10), voxel_size=3.35)

# Initialize with some coherence
for i in range(n_residues):
    atom_idx = i * 4  # CA atom
    voxel_idx = phase_field.get_voxel_indices(positions[atom_idx:atom_idx+1])[0]
    if voxel_idx is not None:
        phase_field.coherence[voxel_idx] = 0.8

# Create analyzer and build graph
analyzer = PatternAnalyzer()
voxel_graph, voxel_map = analyzer._build_voxel_graph(phase_field, positions)

print(f"Positions shape: {positions.shape}")
print(f"Number of voxels: {voxel_graph.number_of_nodes()}")
print(f"Number of edges: {voxel_graph.number_of_edges()}")
print(f"Voxel map: {dict(voxel_map)}")

# Test with torsion states
torsion_states = np.array([2, 2, 3, 2, 3])  # Some helix-like states

# Run full analysis
analysis = analyzer.analyze_template(phase_field, positions, torsion_states)

print(f"\nAnalysis results:")
print(f"  N voxels: {analysis.n_voxels}")
print(f"  N components: {analysis.n_components}")
print(f"  N loops: {analysis.n_loops}")
print(f"  Contact order: {analysis.contact_order:.3f}")
print(f"  Path entropy: {analysis.path_entropy:.3f}")
print(f"  Mobility anisotropy: {analysis.mobility_anisotropy:.3f}")
print(f"  P_ledger: {analysis.p_ledger:.3f}")
print(f"  P_geom: {analysis.p_geom:.3f}") 