#!/usr/bin/env python3
"""Analyze k0 values and their relationship to folding times."""

import numpy as np

# Constants
PHI = 1.618033988749895
TAU_0 = 7.33e-15  # s
kT_room = 0.026  # eV at ~300K

# Test different scenarios
print("k0 Analysis for Protein Folding")
print("="*60)

# For n=20 (Trp-cage)
n = 20
base_k0 = 1 / (8 * TAU_0) * PHI ** (-n / 2)
print(f"\nFor n={n} residues:")
print(f"Base k0: {base_k0:.2e} s^-1")

# Different parameter combinations
scenarios = [
    {"name": "Default", "p_ledger": 0.5, "p_geom": 0.01, "path_eff": 1.0, "mob_factor": 1.0},
    {"name": "Template (old)", "p_ledger": 0.486, "p_geom": 0.618, "path_eff": 1.0, "mob_factor": 1.0},
    {"name": "Template (new)", "p_ledger": 0.486, "p_geom": 1.0, "path_eff": 1.0, "mob_factor": 1.0},
    {"name": "Realistic", "p_ledger": 0.5, "p_geom": 0.1, "path_eff": 0.8, "mob_factor": 0.9},
]

print("\nFolding time predictions:")
print(f"{'Scenario':<20} {'k0 (s^-1)':<12} {'k_fold (s^-1)':<12} {'Mean time (μs)':<15}")
print("-"*60)

for s in scenarios:
    k0 = base_k0 * s["p_ledger"] * s["p_geom"] * s["path_eff"] * s["mob_factor"]
    k_fold = k0 * np.exp(-0.18 / kT_room)  # 0.18 eV barrier
    mean_time_us = 1e6 / k_fold  # Mean of exponential distribution
    
    print(f"{s['name']:<20} {k0:.2e} {k_fold:.2e} {mean_time_us:<15.1f}")

# What k_fold do we need for experimental times?
print("\n\nRequired k_fold for experimental times:")
exp_times = [
    ("Trp-cage", 4.1),
    ("Villin", 0.7),
    ("BBA5", 13.0),
    ("WW domain", 13.0),
]

for name, time_us in exp_times:
    k_fold_needed = 1e6 / time_us
    k0_needed = k_fold_needed / np.exp(-0.18 / kT_room)
    modifier_needed = k0_needed / base_k0
    
    print(f"{name:<15} {time_us:>6.1f} μs → k_fold = {k_fold_needed:.2e} s^-1")
    print(f"{'':15} {'':>9} → k0 = {k0_needed:.2e} s^-1")
    print(f"{'':15} {'':>9} → modifier = {modifier_needed:.2e}")

# Suggest reasonable P_geom values
print("\n\nSuggested P_geom values (assuming p_ledger ≈ 0.5):")
for name, time_us in exp_times:
    k_fold_needed = 1e6 / time_us
    k0_needed = k_fold_needed / np.exp(-0.18 / kT_room)
    modifier_needed = k0_needed / base_k0
    p_geom_suggested = modifier_needed / 0.5  # Assuming p_ledger = 0.5
    
    print(f"{name:<15} P_geom ≈ {p_geom_suggested:.4f}") 