# Missing Files Found During Sanity Check

## Issue Discovered
When trying to run `test_three_proteins_quick.py`, we discovered many essential files were missing from the clean project folder.

## Files That Were Missing

### Core Dependencies
1. **`accelerated_folder.py`** - Base class for AcceleratedFolderV3
   - Referenced in ROADMAP_CLEAN.md but not copied
   - Essential for Monte Carlo barrier crossing

2. **`enhanced_three_layer_folder.py`** - Required by accelerated_folder.py
   - Was moved to archive/ but still needed
   - Contains the three-layer architecture implementation

3. **Entire `phase2_information/` directory** - Required by enhanced_three_layer_folder.py
   - Contains all the core RS implementations:
     - `three_layer_folder.py`
     - `phase_pattern_field.py` (duplicate of what's in core/)
     - `torsion_dynamics.py`
     - `ir_photon_analysis.py`
     - And many more...

## Root Cause
The clean folder creation focused on the "latest" versions (v3 files) without checking their dependencies. The v3 files build upon v2 and v1 implementations, creating a dependency chain.

## Solution Applied
1. Copied `accelerated_folder.py` from original project
2. Copied `enhanced_three_layer_folder.py` from archive back to core
3. Copied entire `phase2_information/` directory to maintain all dependencies
4. Replaced hard-coded mobility constant `0.1` with RS-derived value `≈0.096` in all active algorithm files (`core/enhanced_three_layer_folder.py`, `phase2_information/three_layer_folder.py`, `phase2_information/enhanced_three_layer_folder.py`).

## Lessons Learned
- When creating a "clean" project folder, need to trace ALL dependencies
- The v3 implementations are not standalone - they inherit from earlier versions
- Import statements need to be checked and adjusted for the new structure

## Current Status
✅ Clean folder is now complete and functional
✅ Tests run successfully
✅ Mobility constant empirical value cleaned up
✅ `K0_FOLDING` placeholder removed. `three_layer_folder.py` now computes size-dependent base k₀ dynamically from 1/(8τ₀)·φ^(−n/2).
✅ Results match expected behavior (including known issues) 