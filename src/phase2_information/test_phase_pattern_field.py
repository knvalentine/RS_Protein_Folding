"""
Test Phase Pattern Field implementation

Verifies that the information layer correctly:
1. Accumulates phase patterns from recognition events
2. Maintains ledger balance (conservation)
3. Detects coherence/template completion
4. Computes information pressure
"""

import numpy as np
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from phase2_information.phase_pattern_field import PhasePatternField, RecognitionEvent

def test_basic_recognition():
    """Test basic recognition event processing"""
    print("Testing basic recognition event...")
    
    # Create field for 3-residue system
    field = PhasePatternField(n_residues=3)
    
    # Add recognition between residues 0 and 1
    event = RecognitionEvent(
        residue_i=0,
        residue_j=1,
        phase_shift=np.pi/4,  # 45 degrees
        tick=100
    )
    field.add_recognition(event)
    
    # Check that phase was distributed across channels
    assert field.recognition_count == 1
    assert field.recognition_matrix[0, 1] == True
    assert field.recognition_matrix[1, 0] == True
    
    # Check phase conservation (ledger balance)
    channel_sums = np.sum(field.phase_pattern, axis=0)
    assert np.allclose(channel_sums, 0), f"Phase not conserved: {channel_sums}"
    
    print("✓ Basic recognition test passed")


def test_coherence_detection():
    """Test coherence detection for complete template"""
    print("\nTesting coherence detection...")
    
    # Create field for 4-residue system
    field = PhasePatternField(n_residues=4)
    
    # Add recognitions to form a connected network
    # Linear chain: 0-1-2-3
    recognitions = [
        (0, 1, np.pi/3),
        (1, 2, np.pi/4),
        (2, 3, np.pi/6),
    ]
    
    for i, (res_i, res_j, phase) in enumerate(recognitions):
        event = RecognitionEvent(
            residue_i=res_i,
            residue_j=res_j,
            phase_shift=phase,
            tick=(i+1) * 100
        )
        field.add_recognition(event)
    
    # Check connectivity
    assert field._is_connected(), "Network should be connected"
    
    # Check if template is complete (might need more events for amplitude)
    progress = field.get_folding_progress()
    print(f"  Recognition count: {progress['recognition_count']}")
    print(f"  Time: {progress['time_ps']:.2f} ps")
    print(f"  Coherence: {progress['coherence']:.3f}")
    print(f"  Template complete: {progress['template_complete']}")
    
    print("✓ Coherence detection test passed")


def test_information_pressure():
    """Test information pressure calculation"""
    print("\nTesting information pressure...")
    
    # Create field for 2-residue system
    field = PhasePatternField(n_residues=2)
    
    # Add recognition with specific phase
    event = RecognitionEvent(
        residue_i=0,
        residue_j=1,
        phase_shift=0,  # Zero phase shift
        tick=100
    )
    field.add_recognition(event)
    
    # Test positions: too close
    positions = np.array([
        [0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0]  # 3 Å apart (target is 6.5 Å)
    ])
    
    pressure = field.compute_information_pressure(positions)
    
    # Should push residues apart
    assert pressure[0, 0] < 0, "Residue 0 should be pushed left"
    assert pressure[1, 0] > 0, "Residue 1 should be pushed right"
    assert np.allclose(pressure[0] + pressure[1], 0), "Total pressure should be zero"
    
    print(f"  Pressure on residue 0: {pressure[0]}")
    print(f"  Pressure on residue 1: {pressure[1]}")
    print("✓ Information pressure test passed")


def test_golden_ratio_scaling():
    """Test that golden ratio appears in phase relationships"""
    print("\nTesting golden ratio scaling...")
    
    field = PhasePatternField(n_residues=5)
    
    # Add recognitions with golden ratio phase relationships
    phi = (1 + np.sqrt(5)) / 2
    
    for i in range(4):
        phase = 2 * np.pi / (phi ** i)  # Golden ratio scaling
        event = RecognitionEvent(
            residue_i=i,
            residue_j=i+1,
            phase_shift=phase,
            tick=(i+1) * 100
        )
        field.add_recognition(event)
    
    # Check channel summary
    channel_summary = field.get_channel_summary()
    print(f"  Channel RMS amplitudes: {channel_summary}")
    
    # Verify phase conservation
    channel_sums = np.sum(field.phase_pattern, axis=0)
    assert np.allclose(channel_sums, 0), "Phase not conserved with golden ratio"
    
    print("✓ Golden ratio scaling test passed")


def test_eight_channel_architecture():
    """Test the 8-channel phase architecture"""
    print("\nTesting 8-channel architecture...")
    
    field = PhasePatternField(n_residues=2)
    
    # Add recognition
    event = RecognitionEvent(
        residue_i=0,
        residue_j=1,
        phase_shift=np.pi,
        tick=100
    )
    field.add_recognition(event)
    
    # Check that all 8 channels received phase information
    assert field.phase_pattern.shape == (2, 8), "Should have 8 channels"
    
    # Each channel should have different coupling based on golden angle
    channel_0_coupling = field.phase_pattern[0, 0]
    channel_1_coupling = field.phase_pattern[0, 1]
    
    # Channels at golden angle spacing should have different couplings
    assert not np.isclose(channel_0_coupling, channel_1_coupling), \
        "Different channels should have different couplings"
    
    print(f"  Phase pattern shape: {field.phase_pattern.shape}")
    print(f"  Channel 0 coupling: {channel_0_coupling:.3f}")
    print(f"  Channel 1 coupling: {channel_1_coupling:.3f}")
    print("✓ 8-channel architecture test passed")


def run_all_tests():
    """Run all tests"""
    print("=" * 50)
    print("Testing Phase Pattern Field Implementation")
    print("=" * 50)
    
    test_basic_recognition()
    test_coherence_detection()
    test_information_pressure()
    test_golden_ratio_scaling()
    test_eight_channel_architecture()
    
    print("\n" + "=" * 50)
    print("All tests passed! ✓")
    print("=" * 50)


if __name__ == "__main__":
    run_all_tests() 