"""
make_all_figures.py - Reproduce all figures from the paper
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('..')
from sct_core import CAR_predictions


def figure_1_car_predictions():
    """Figure 1: CAR predictions vs. ΛCDM"""
    preds = CAR_predictions()
    
    # Create bar chart comparing predictions
    categories = ['r_d (Mpc)', 'H₀ (km/s/Mpc)', 'S₈', 'IA bias']
    car_values = [preds['r_d_Mpc'], preds['H0_km_s_Mpc'], preds['S8'], preds['IA_bias']]
    lcdm_values = [150.0, 67.4, 0.832, 1.0]
    
    x = np.arange(len(categories))
    width = 0.35
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(x - width/2, car_values, width, label='CAR', color='blue')
    ax.bar(x + width/2, lcdm_values, width, label='ΛCDM', color='gray', alpha=0.7)
    
    ax.set_ylabel('Value')
    ax.set_title('CAR Predictions vs. ΛCDM')
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figure_1_car_predictions.png', dpi=150)
    plt.close()


def figure_2_residuals():
    """Figure 2: Residual plots for each dataset"""
    datasets = ['DESI-DR2 BAO', 'Planck PR4', 'DES-Y6', 'HSC+KiDS']
    chi2_dof = [0.931, 1.028, 0.982, 0.994]
    colors = ['blue', 'green', 'red', 'purple']
    
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(datasets, chi2_dof, color=colors, alpha=0.7)
    ax.axhline(y=1.0, color='black', linestyle='--', label='ΛCDM baseline (χ²/dof=1)')
    ax.set_ylabel('χ² / dof')
    ax.set_title('Goodness of Fit: CAR Predictions')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figure_2_residuals.png', dpi=150)
    plt.close()


def figure_3_s8_comparison():
    """Figure 3: S₈ comparison across surveys"""
    surveys = ['DES-Y6', 'HSC-Y3', 'KiDS-DR5', 'Planck', 'CAR']
    s8_values = [0.780, 0.776, 0.788, 0.832, 0.7838]   # v4.8.1: KiDS 0.815→0.788, S8 SCT 0.783→0.7838
    s8_errors = [0.012, 0.032, 0.016, 0.013, 0.015]
    colors = ['red', 'orange', 'green', 'blue', 'black']
    
    fig, ax = plt.subplots(figsize=(10, 6))
    for i, (survey, val, err, color) in enumerate(zip(surveys, s8_values, s8_errors, colors)):
        ax.errorbar(i, val, yerr=err, fmt='o', color=color, capsize=5, markersize=10)
        ax.text(i, val + err + 0.005, survey, ha='center', fontsize=9)
    
    ax.axhline(y=0.783, color='black', linestyle='--', label='CAR Prediction (0.783)')
    ax.set_ylabel('S₈')
    ax.set_xticks([])
    ax.set_title('S₈ Measurements vs. CAR Prediction')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figure_3_s8_comparison.png', dpi=150)
    plt.close()


def figure_4_h0_comparison():
    """Figure 4: H₀ comparison across measurements"""
    measurements = ['SH0ES', 'Planck', 'DESI+Planck', 'CAR']
    h0_values = [73.0, 67.4, 70.4, 70.4]
    h0_errors = [1.0, 0.5, 0.5, 0.5]
    colors = ['red', 'blue', 'green', 'black']
    
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(measurements))
    ax.errorbar(x, h0_values, yerr=h0_errors, fmt='o', color='black', capsize=5, markersize=10)
    
    for i, (name, val, err, color) in enumerate(zip(measurements, h0_values, h0_errors, colors)):
        ax.plot(i, val, 'o', color=color, markersize=12)
        ax.text(i, val + err + 0.5, name, ha='center', fontsize=10)
    
    ax.axhline(y=70.4, color='black', linestyle='--', label='CAR Prediction (70.4)')
    ax.set_ylabel('H₀ (km/s/Mpc)')
    ax.set_xticks([])
    ax.set_title('H₀ Measurements vs. CAR Prediction')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figure_4_h0_comparison.png', dpi=150)
    plt.close()


def figure_5_bayesian_evidence():
    """Figure 5: Bayesian evidence comparison"""
    models = ['ΛCDM', 'w₀w_a', 'EDE', 'CAR']
    delta_lnB = [0, -2.2, -1.7, -3.8]
    colors = ['gray', 'orange', 'green', 'blue']
    
    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(models, delta_lnB, color=colors, alpha=0.7)
    ax.axhline(y=-3.8, color='blue', linestyle='--', linewidth=2, label='CAR Δln B = -3.8')
    ax.set_ylabel('Δln B (relative to ΛCDM)')
    ax.set_title('Bayesian Evidence Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Add text labels
    for bar, val in zip(bars, delta_lnB):
        if val < 0:
            ax.text(bar.get_x() + bar.get_width()/2, val - 0.2, f'{val:.1f}', 
                    ha='center', va='top', fontsize=10)
    
    plt.tight_layout()
    plt.savefig('figure_5_bayesian_evidence.png', dpi=150)
    plt.close()


def figure_6_kill_switch():
    """Figure 6: Euclid Y1 kill-switch visualization"""
    # Simulated Euclid Y1 measurement scenarios
    car_prediction = 0.781
    car_error = 0.009
    kill_threshold = 0.81
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # CAR prediction region
    car_region = np.array([car_prediction - car_error, car_prediction + car_error])
    ax.axvspan(car_region[0], car_region[1], alpha=0.3, color='green', label='CAR 1σ region')
    
    # Kill region
    ax.axvspan(kill_threshold, 0.84, alpha=0.3, color='red', label='Kill zone (S₈ > 0.81)')
    
    # Vertical line at CAR prediction
    ax.axvline(x=car_prediction, color='green', linestyle='-', linewidth=2, label='CAR Prediction')
    
    # Kill threshold line
    ax.axvline(x=kill_threshold, color='red', linestyle='--', linewidth=2, label='Kill threshold (3σ)')
    
    ax.set_xlabel('S₈')
    ax.set_ylabel('Probability Density')
    ax.set_title('Euclid Y1: CAR Kill-Switch')
    ax.legend()
    ax.set_xlim(0.75, 0.84)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figure_6_kill_switch.png', dpi=150)
    plt.close()


if __name__ == "__main__":
    print("Generating figures...")
    figure_1_car_predictions()
    figure_2_residuals()
    figure_3_s8_comparison()
    figure_4_h0_comparison()
    figure_5_bayesian_evidence()
    figure_6_kill_switch()
    print("Done. Figures saved to current directory.")