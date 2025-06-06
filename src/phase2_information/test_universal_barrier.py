"""RS-truth unit test: every analyzer must report ΔC = 0.18 eV (2 coins)."""

import math
from pattern_analyzer import PatternAnalyzer, E_COH as E_COH_BASE
from pattern_analyzer_v2 import PatternAnalyzerV2, E_COH as E_COH_V2


def test_pattern_analyzer_barrier():
    pa = PatternAnalyzer()
    barrier = pa._calculate_barrier(set([0, 1, 2]), n_components=3)  # even with complex input
    assert barrier == 2, "PatternAnalyzer must return 2 coins universally"


def test_pattern_analyzer_v2_barrier():
    pav2 = PatternAnalyzerV2()
    barrier_ev = pav2._calculate_effective_barrier(
        coherence_length=1.0,
        info_flow=0.1,
        rec_density=100.0,
        frustration=0.9,
    )
    assert math.isclose(barrier_ev, 2 * E_COH_V2, rel_tol=1e-6), "Barrier must be fixed at 0.18 eV"


if __name__ == "__main__":
    test_pattern_analyzer_barrier()
    test_pattern_analyzer_v2_barrier()
    print("All universal barrier tests passed.") 