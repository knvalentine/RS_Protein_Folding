# Remaining ROADMAP_CLEAN.md Fixes Needed

## Line 221 needs to be changed from:
```
  - [ ] Derive barrier and prefactors from topology
```

## To:
```
  - [ ] Derive prefactors from topology (barrier is always 0.18 eV)
```

## Why:
The barrier is axiomatically fixed at 2 E_coh = 0.18 eV by Recognition Science principles. Only the prefactor k₀ varies with protein topology.

## Already Fixed:
✅ Phase 2.7 now correctly states barrier is always 0.18 eV
✅ Phase 2.8 bullets correctly say "barrier fixed at 2 coins"
✅ Pattern analyzers now return universal barrier
✅ Tests verify barrier is always 2 coins

## Summary:
The roadmap is 99% correct. Just need to fix line 221 to not mention deriving the barrier. 