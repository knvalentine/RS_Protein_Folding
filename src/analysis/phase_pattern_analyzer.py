"""
Phase Pattern Analyzer - Extract predictive features from template phase patterns

This module implements advanced analysis of the PhasePatternField to find
correlations between template-stage patterns and final β-sheet formation.

Key idea: The 8-channel phase architecture may encode future structural
information that we're not currently extracting.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'core'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'phase2_information'))
from phase_pattern_field import PhasePatternField

class PhasePatternAnalyzer:
    """Advanced analysis of phase patterns for β-sheet prediction."""
    
    def __init__(self, phase_field: PhasePatternField):
        self.phase_field = phase_field
        self.n_residues = phase_field.n_residues
        
    def analyze_phase_correlations(self) -> Dict[str, float]:
        """
        Analyze phase correlations that might predict β-sheet formation.
        
        Returns:
            Dictionary of correlation metrics
        """
        results = {}
        
        # 1. Channel-wise phase coherence
        channel_coherence = []
        for channel in range(8):
            phases = self.phase_field.phase_pattern[:, channel]
            if np.any(phases):
                # Kuramoto order parameter
                z = np.mean(np.exp(1j * phases))
                channel_coherence.append(np.abs(z))
            else:
                channel_coherence.append(0.0)
        
        results['channel_coherence'] = channel_coherence
        results['max_coherence'] = max(channel_coherence)
        results['coherence_variance'] = np.var(channel_coherence)
        
        # 2. Phase gradient analysis (might indicate strand directions)
        phase_gradients = []
        for i in range(self.n_residues - 1):
            for channel in range(8):
                grad = self.phase_field.phase_pattern[i+1, channel] - \
                       self.phase_field.phase_pattern[i, channel]
                phase_gradients.append(grad)
        
        results['mean_gradient'] = np.mean(np.abs(phase_gradients))
        results['gradient_variance'] = np.var(phase_gradients)
        
        # 3. Long-range phase correlations (β-sheets have long-range order)
        long_range_corr = self._compute_long_range_correlations()
        results['long_range_correlation'] = long_range_corr
        
        # 4. Phase frustration index (competing phase relationships)
        frustration = self._compute_phase_frustration()
        results['phase_frustration'] = frustration
        
        # 5. Channel coupling strength
        coupling = self._compute_channel_coupling()
        results['channel_coupling'] = coupling
        
        return results
    
    def _compute_long_range_correlations(self) -> float:
        """
        Compute correlations between distant residues in phase space.
        β-sheets should show stronger long-range correlations.
        """
        correlations = []
        
        for separation in range(5, min(20, self.n_residues // 2)):
            for i in range(self.n_residues - separation):
                j = i + separation
                
                # Compute phase correlation across all channels
                phase_i = self.phase_field.phase_pattern[i, :]
                phase_j = self.phase_field.phase_pattern[j, :]
                
                if np.any(phase_i) and np.any(phase_j):
                    # Complex correlation
                    z_i = np.mean(np.exp(1j * phase_i))
                    z_j = np.mean(np.exp(1j * phase_j))
                    corr = np.real(z_i * np.conj(z_j))
                    correlations.append(corr)
        
        return np.mean(correlations) if correlations else 0.0
    
    def _compute_phase_frustration(self) -> float:
        """
        Compute phase frustration - competing phase relationships
        that might indicate β-sheet propensity.
        """
        frustration = 0.0
        count = 0
        
        # Check triangles in recognition network
        recog_matrix = self.phase_field.recognition_matrix
        
        for i in range(self.n_residues):
            for j in range(i+1, self.n_residues):
                for k in range(j+1, self.n_residues):
                    if recog_matrix[i,j] and recog_matrix[j,k] and recog_matrix[i,k]:
                        # Triangle exists - check phase consistency
                        edge_ij = self.phase_field.voxel_edges.get((i,j), 0)
                        edge_jk = self.phase_field.voxel_edges.get((j,k), 0)
                        edge_ik = self.phase_field.voxel_edges.get((i,k), 0)
                        
                        # Phase around loop should sum to 2π*n
                        loop_phase = edge_ij + edge_jk - edge_ik
                        deviation = abs(loop_phase % (2*np.pi))
                        frustration += min(deviation, 2*np.pi - deviation)
                        count += 1
        
        return frustration / count if count > 0 else 0.0
    
    def _compute_channel_coupling(self) -> float:
        """
        Compute coupling between different phase channels.
        Strong coupling might indicate structural constraints.
        """
        coupling_matrix = np.zeros((8, 8))
        
        for ch1 in range(8):
            for ch2 in range(ch1+1, 8):
                phases1 = self.phase_field.phase_pattern[:, ch1]
                phases2 = self.phase_field.phase_pattern[:, ch2]
                
                if np.any(phases1) and np.any(phases2):
                    # Compute phase locking
                    phase_diff = phases1 - phases2
                    z = np.mean(np.exp(1j * phase_diff))
                    coupling_matrix[ch1, ch2] = np.abs(z)
        
        return np.mean(coupling_matrix[coupling_matrix > 0])
    
    def predict_beta_registries(self) -> int:
        """
        Predict number of β-registries based on phase pattern analysis.
        
        This is our key innovation - using information theory to predict
        future structural features from current phase patterns.
        """
        metrics = self.analyze_phase_correlations()
        
        # Empirical prediction based on phase metrics
        # TODO: This needs calibration with known proteins
        score = 0.0
        
        # Long-range correlations indicate β-sheet potential
        score += metrics['long_range_correlation'] * 10
        
        # High frustration suggests competing structures (often β)
        score += metrics['phase_frustration'] * 5
        
        # Low coherence variance suggests ordered structure
        score += (1 - metrics['coherence_variance']) * 3
        
        # Channel coupling indicates structural constraints
        score += metrics['channel_coupling'] * 2
        
        # Threshold-based prediction
        if score > 8.0:
            return 2  # Likely multiple β-sheets
        elif score > 4.0:
            return 1  # Likely single β-sheet
        else:
            return 0  # Likely all-α
    
    def get_predictive_features(self) -> Dict[str, float]:
        """
        Extract all features that might be predictive of final structure.
        
        These features can be used for machine learning or empirical
        correlation with experimental folding times.
        """
        features = self.analyze_phase_correlations()
        
        # Add eigenvalue analysis of phase correlation matrix
        phase_corr_matrix = np.corrcoef(self.phase_field.phase_pattern.T)
        if not np.any(np.isnan(phase_corr_matrix)):
            eigenvalues = np.linalg.eigvalsh(phase_corr_matrix)
            features['largest_eigenvalue'] = eigenvalues[-1]
            features['eigenvalue_gap'] = eigenvalues[-1] - eigenvalues[-2]
            features['eigenvalue_entropy'] = -np.sum(eigenvalues * np.log(eigenvalues + 1e-10))
        
        # Add recognition network metrics
        recog_density = np.sum(self.phase_field.recognition_matrix) / (self.n_residues**2)
        features['recognition_density'] = recog_density
        
        # Predict β-registries
        features['predicted_beta_registries'] = self.predict_beta_registries()
        
        return features 