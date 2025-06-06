#!/usr/bin/env python3
"""Test RS-derived estimators with various graph topologies."""

import numpy as np
import networkx as nx
from pattern_analyzer import PatternAnalyzer

# Test different graph topologies
test_graphs = {
    "linear_chain": nx.path_graph(5),
    "cycle": nx.cycle_graph(6),
    "grid_2x3": nx.grid_2d_graph(2, 3),
    "complete_k4": nx.complete_graph(4),
    "barbell": nx.barbell_graph(3, 2),
}

analyzer = PatternAnalyzer()

print("RS-Derived Estimator Tests")
print("="*60)

for name, G in test_graphs.items():
    # Convert to 3D coordinates for consistency
    G_3d = nx.Graph()
    for i, node in enumerate(G.nodes()):
        if isinstance(node, tuple):
            # Already has coordinates (like grid)
            G_3d.add_node((node[0], node[1], 0))
        else:
            # Linear layout
            G_3d.add_node((i, 0, 0))
    
    # Copy edges
    for u, v in G.edges():
        if isinstance(u, tuple):
            u_3d = (u[0], u[1], 0)
            v_3d = (v[0], v[1], 0)
        else:
            u_idx = list(G.nodes()).index(u)
            v_idx = list(G.nodes()).index(v)
            u_3d = (u_idx, 0, 0)
            v_3d = (v_idx, 0, 0)
        G_3d.add_edge(u_3d, v_3d)
    
    print(f"\n{name.upper()}")
    print(f"  Nodes: {G_3d.number_of_nodes()}")
    print(f"  Edges: {G_3d.number_of_edges()}")
    print(f"  Components: {nx.number_connected_components(G_3d)}")
    
    # Calculate metrics
    path_entropy = analyzer._estimate_path_entropy(G_3d)
    mobility = analyzer._estimate_mobility_anisotropy(G_3d)
    n_loops = analyzer._count_independent_loops(G_3d)
    
    print(f"  Loops: {n_loops}")
    print(f"  Path entropy: {path_entropy:.4f}")
    print(f"  Mobility anisotropy: {mobility:.4f}")
    
    # Calculate what k0 modifier would be
    p_ledger = 1.618033988749895 ** (-1 / 2)  # Single component
    p_geom = 1.618033988749895 ** (-n_loops / 2)
    path_efficiency = 1.0 - path_entropy
    mobility_factor = 1.0 / (1.0 + mobility)
    
    k0_modifier = p_ledger * p_geom * path_efficiency * mobility_factor
    print(f"  k0 modifier: {k0_modifier:.6f}")

# Test the full formula with realistic values
print("\n" + "="*60)
print("FULL k0 CALCULATION TEST")
print("="*60)

PHI = 1.618033988749895
TAU_0 = 7.33e-15
n_residues = 20

# Base k0
base_k0 = 1 / (8 * TAU_0) * PHI ** (-n_residues / 2)
print(f"Base k0 (n=20): {base_k0:.2e} s^-1")

# With default parameters
default_k0 = base_k0 * 0.5 * 0.01  # p_ledger=0.5, p_geom=0.01
print(f"Default k0: {default_k0:.2e} s^-1")

# With template parameters (example from cycle graph)
p_ledger = PHI ** (-1/2)  # 1 component
p_geom = PHI ** (-1/2)    # 1 loop
path_entropy = 0.5        # Some redundancy
mobility = 0.0            # Isotropic

path_efficiency = 1.0 - path_entropy
mobility_factor = 1.0 / (1.0 + mobility)
template_k0 = base_k0 * p_ledger * p_geom * path_efficiency * mobility_factor

print(f"\nTemplate-driven k0:")
print(f"  p_ledger: {p_ledger:.4f}")
print(f"  p_geom: {p_geom:.4f}")
print(f"  path_efficiency: {path_efficiency:.4f}")
print(f"  mobility_factor: {mobility_factor:.4f}")
print(f"  Final k0: {template_k0:.2e} s^-1")
print(f"  Ratio to default: {template_k0/default_k0:.2f}x") 