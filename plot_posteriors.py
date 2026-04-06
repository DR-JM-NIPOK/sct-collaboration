#!/usr/bin/env python
"""
plot_posteriors.py - Generate posterior plots from chains
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os


def plot_posteriors(chain_file, output_file):
    """
    Generate posterior plots from PolyChord chains.
    
    Parameters
    ----------
    chain_file : str
        Path to PolyChord chain file
    output_file : str
        Path to output plot file
    """
    try:
        # Read chains
        chains = np.loadtxt(chain_file)
        
        # Extract parameters (assuming standard order)
        # Columns: Omega_m, R_b, logL, weight, etc.
        Omega_m = chains[:, 0]
        R_b = chains[:, 1]
        logL = chains[:, -1]
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Omega_m histogram
        axes[0, 0].hist(Omega_m, bins=50, density=True, alpha=0.7, color='blue')
        axes[0, 0].axvline(np.mean(Omega_m), color='red', linestyle='--', label=f'Mean: {np.mean(Omega_m):.3f}')
        axes[0, 0].axvline(np.percentile(Omega_m, 16), color='gray', linestyle=':')
        axes[0, 0].axvline(np.percentile(Omega_m, 84), color='gray', linestyle=':')
        axes[0, 0].set_xlabel('Ω_m')
        axes[0, 0].set_ylabel('Density')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # R_b histogram
        axes[0, 1].hist(R_b, bins=50, density=True, alpha=0.7, color='green')
        axes[0, 1].axvline(np.mean(R_b), color='red', linestyle='--', label=f'Mean: {np.mean(R_b):.3f}')
        axes[0, 1].axvline(0.2545, color='black', linestyle='-', label='Derived: 0.2545 (Paper 17 v4.8 Section 11.6)')
        axes[0, 1].set_xlabel('R_b')
        axes[0, 1].set_ylabel('Density')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # Omega_m vs R_b scatter
        axes[1, 0].scatter(Omega_m, R_b, alpha=0.1, s=1, color='purple')
        axes[1, 0].set_xlabel('Ω_m')
        axes[1, 0].set_ylabel('R_b')
        axes[1, 0].grid(True, alpha=0.3)
        
        # LogL vs Omega_m
        axes[1, 1].scatter(Omega_m, logL, alpha=0.1, s=1, color='orange')
        axes[1, 1].set_xlabel('Ω_m')
        axes[1, 1].set_ylabel('log Likelihood')
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.suptitle('CAR Model Posteriors')
        plt.tight_layout()
        plt.savefig(output_file, dpi=150)
        plt.close()
        
        print(f"Posterior plot saved to {output_file}")
        
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Input chain file')
    parser.add_argument('--output', required=True, help='Output plot file')
    args = parser.parse_args()
    
    plot_posteriors(args.input, args.output)