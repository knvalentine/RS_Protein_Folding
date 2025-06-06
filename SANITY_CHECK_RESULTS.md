# Sanity Check Results

## Test Performed
Attempted to run `test_three_proteins_quick.py` with only the files in the clean folder.

## Issues Found

### 1. Missing Core Dependencies
- `accelerated_folder.py` was not copied (but referenced in docs)
- `enhanced_three_layer_folder.py` was in archive/ but still needed
- Entire `phase2_information/` directory was missing

### 2. Import Path Issues
- Tests assumed flat import structure
- Had to add sys.path modifications to find modules

### 3. Documentation vs Reality Mismatch
- ROADMAP_CLEAN.md references `accelerated_folder.py` multiple times
- IMPORTANT_FILES.md lists it as 3-star importance
- But it wasn't in the clean folder initially

## Resolution
1. Copied missing files from original project
2. Fixed import paths
3. Replaced empirical mobility constant `0.1` with RS-derived `≈0.096`
4. Test now runs successfully

## Test Results
After fixing dependencies:
- Trp-cage: 0.0 μs (357× too fast) - β-registry detection issue
- BBA5: 0.0 μs (359× too fast) - β-registry detection issue  
- Villin: 0.4 μs vs 0.7 μs experimental ✅ Within factor of 2!

This confirms the β-registry detection problem documented in ROADMAP_CLEAN.md.

## Recommendations
1. **For future clean folders**: Use dependency analysis tools or trace imports
2. **For this project**: The folder now contains all needed files
3. **For documentation**: Update to clarify V3 depends on V2/V1 implementations
4. Placeholder `K0_FOLDING` removed—`three_layer_folder.py` now derives base k₀ dynamically. 