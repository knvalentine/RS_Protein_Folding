#!/usr/bin/env python3
"""Debug the pattern analyzer to see what values it's producing."""

import numpy as np
import networkx as nx
from pattern_analyzer import PatternAnalyzer

# Create a simple test graph
G = nx.Graph()
# Add a simple 3x3 grid of voxels
for i in range(3):
    for j in range(3):
        G.add_node((i, j, 0))

# Connect adjacent voxels
for i in range(3):
    for j in range(3):
        node = (i, j, 0)
        # Right neighbor
        if i < 2:
            G.add_edge(node, (i+1, j, 0))
        # Down neighbor  
        if j < 2:
            G.add_edge(node, (i, j+1, 0))

print("Test graph:")
print(f"Nodes: {G.number_of_nodes()}")
print(f"Edges: {G.number_of_edges()}")
print(f"Connected components: {nx.number_connected_components(G)}")

# Test the estimator functions
analyzer = PatternAnalyzer()

# Test path entropy
path_entropy = analyzer._estimate_path_entropy(G)
print(f"\nPath entropy: {path_entropy}")

# Test mobility anisotropy
mobility = analyzer._estimate_mobility_anisotropy(G)
print(f"Mobility anisotropy: {mobility}")

# Test with a linear chain (should have higher anisotropy)
G2 = nx.path_graph(5)
# Convert to 3D coordinates
G2_3d = nx.Graph()
for i in range(5):
    G2_3d.add_node((i, 0, 0))
for i in range(4):
    G2_3d.add_edge((i, 0, 0), (i+1, 0, 0))

print("\n\nLinear chain graph:")
path_entropy2 = analyzer._estimate_path_entropy(G2_3d)
mobility2 = analyzer._estimate_mobility_anisotropy(G2_3d)
print(f"Path entropy: {path_entropy2}")
print(f"Mobility anisotropy: {mobility2}")

# Test the k0 calculation
PHI = 1.618033988749895
TAU_0 = 7.33e-15
n_residues = 20

# Default values
p_ledger_default = 0.5
p_geom_default = 0.01
base_k0 = 1 / (8 * TAU_0) * PHI ** (-n_residues / 2)
k0_default = base_k0 * p_ledger_default * p_geom_default
print(f"\n\nDefault k0 calculation (n=20):")
print(f"Base k0: {base_k0:.2e}")
print(f"With defaults: {k0_default:.2e}")

# With template values (example)
p_ledger = 0.486
p_geom = 0.618
k0_template = base_k0 * p_ledger * p_geom * path_entropy * (1 / (1 + mobility))
print(f"\nTemplate-driven k0:")
print(f"p_ledger: {p_ledger}")
print(f"p_geom: {p_geom}")
print(f"path_entropy: {path_entropy}")
print(f"1/(1+mobility): {1/(1+mobility)}")
print(f"Final k0: {k0_template:.2e}")

# Check if k0 is zero
if k0_template == 0:
    print("\nWARNING: k0 is zero! This will cause divide-by-zero error.")
    print("Path entropy might be returning 0 for well-connected graphs.") 