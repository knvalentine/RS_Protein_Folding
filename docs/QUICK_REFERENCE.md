# Recognition Science Quick Reference

## Fundamental Constants
- **E_coh** = 0.090 eV (coherence energy)
- **τ₀** = 58.6 fs (eight-beat cycle time)
- **φ** = 1.618... (golden ratio)
- **Voxel size** = 3.35 Å (2φ Bohr radii)

## Key Formulas

### Folding Rate
```
k_fold = k₀ · exp(-ΔG/kT)
```

### Universal Barrier
```
ΔG = 2 · E_coh = 0.18 eV (ALL proteins)
```

### Prefactor
```
k₀ = (1/8τ₀) · φ^(-n/2) · P_ledger · P_geom · (1-S_path) · mobility_factor
```

### Geometric Factor
```
P_geom = φ^(-loops/2) · φ^(-√φ·CO) · φ^(L/2) · φ^(-R)
```
Where:
- loops = independent cycles in voxel graph
- CO = contact order (normalized)
- L = number of helix axes (≥6 residues)
- R = number of β-registries

### Path Entropy
```
S_path = (1/N_pairs) · Σ_{u<v} ln[N_sp(u,v)] / ln(φ)
```
- N_sp(u,v) = number of shortest paths between voxels u,v
- Normalized to [0,1] range

### Mobility Anisotropy
```
Anisotropy = (λ_max/λ_min) - 1
```
- From eigenvalues of voxel cloud covariance

## Two Timescales
1. **Template Formation**: 5-20 ps
   - Phase patterns organize
   - Backbone connectivity established
   - Information layer complete

2. **Physical Folding**: 0.1-100 μs
   - Barrier crossing (0.18 eV)
   - Secondary structure formation
   - Final native state

## Torsion Glyphs (9 states)
```
Glyph 0: φ=-120°, ψ=-120° (sheet)
Glyph 1: φ=-120°, ψ=  0°
Glyph 2: φ=-120°, ψ=+120° (helix)
Glyph 3: φ=  0°, ψ=-120° (helix)
Glyph 4: φ=  0°, ψ=  0°
Glyph 5: φ=  0°, ψ=+120°
Glyph 6: φ=+120°, ψ=-120°
Glyph 7: φ=+120°, ψ=  0°
Glyph 8: φ=+120°, ψ=+120° (sheet)
```

## Recognition Rules
1. **Discrete Events**: No continuous forces
2. **Phase Modulation**: P_recognize = |cos(Δφ/2)|²
3. **Conservation**: Coins, momentum, energy preserved
4. **Eight-Beat Cycle**: Events occur on τ₀ ticks

## Current Limitations
- Template analysis at ps timescale misses final β-sheets
- L and R often = 0 at template stage
- Need better topology prediction

## File Organization
```
src/core/
  recognition_dynamics_v4.py    # Core RS mechanics
  phase_pattern_field.py        # Information layer
  pattern_analyzer.py           # Template analysis
  accelerated_folder_v3.py      # MC acceleration

theory/
  Deeper Understanding.txt      # Key derivations
  Protein-Full-2.tex           # Complete theory
```

## Running Tests
```bash
cd src/tests
python test_three_proteins_quick.py  # Quick validation
python test_ten_protein_suite.py     # Full suite
```

## Debug Flags
- Set `debug=True` in pattern_analyzer for detailed output
- Check recognition matrix with `phase_field.recognition_matrix`
- Monitor template completion with `phase_field.get_template_completeness()` 