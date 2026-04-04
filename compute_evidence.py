#!/usr/bin/env python
"""
compute_evidence.py - Compute Bayesian evidence from PolyChord chains
"""

import numpy as np
import sys
import os


def compute_evidence(chain_file, output_file):
    """
    Compute Bayesian evidence from PolyChord output.
    
    Parameters
    ----------
    chain_file : str
        Path to PolyChord chain file (e.g., chains/car_combined.txt)
    output_file : str
        Path to output evidence file
    """
    try:
        # Read the evidence from PolyChord's .stats file
        stats_file = chain_file.replace('.txt', '.stats')
        
        if os.path.exists(stats_file):
            with open(stats_file, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if 'logZ' in line:
                        # Parse log evidence
                        parts = line.split()
                        logZ = float(parts[1])
                        logZ_err = float(parts[2])
                        
                        with open(output_file, 'w') as out:
                            out.write(f"log_evidence: {logZ:.4f} ± {logZ_err:.4f}\n")
                            out.write(f"evidence: {np.exp(logZ):.2e}\n")
                        return
        else:
            # Fallback: compute from chain
            chains = np.loadtxt(chain_file)
            # Simple harmonic mean estimator (not rigorous, for illustration)
            logL = chains[:, -1]  # Last column is log-likelihood
            logZ_est = np.log(np.mean(np.exp(logL - np.max(logL)))) + np.max(logL)
            
            with open(output_file, 'w') as out:
                out.write(f"log_evidence (estimated): {logZ_est:.4f}\n")
                out.write("Note: Use .stats file for accurate evidence\n")
                
    except Exception as e:
        print(f"Error: {e}")
        with open(output_file, 'w') as out:
            out.write(f"ERROR: {e}\n")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Input chain file')
    parser.add_argument('--output', required=True, help='Output evidence file')
    args = parser.parse_args()
    
    compute_evidence(args.input, args.output)