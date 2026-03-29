"""
chains/plot_posteriors.py
==========================
Generate publication-quality posterior contour plots from PolyChord
chain output using GetDist.

Produces:
    posteriors_car_combined.pdf    — CAR Ω_m, R_b, H₀, S₈ posteriors
    posteriors_comparison.pdf      — CAR vs ΛCDM S₈–Ω_m contours
    tension_plot.pdf               — 1D marginals for H₀, S₈, r_d

Author : DR JM NIPOK | License: GPL-3.0
"""

import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from sct_core import CAR_predictions


def plot_car_posteriors(chain_dir: str, output_path: str,
                        model: str = 'car') -> None:
    """
    Generate posterior triangle plot from PolyChord chains.

    Parameters
    ----------
    chain_dir  : directory containing PolyChord .txt chains
    output_path: output PDF path
    model      : 'car' or 'lcdm'
    """
    try:
        from getdist import plots, MCSamples
        import getdist.mcsamples as mcs
    except ImportError:
        print("  [SKIP] getdist not installed. pip install getdist")
        _plot_fallback(output_path, model)
        return

    # Find chain file
    chain_file = None
    for root, dirs, files in os.walk(chain_dir):
        for f in files:
            if f.endswith('.txt') and model in f:
                chain_file = os.path.join(root, f)
                break

    if chain_file is None or not os.path.exists(chain_file):
        print(f"  [INFO] No chain file found in {chain_dir}; generating demo plot.")
        _plot_fallback(output_path, model)
        return

    # Load chains
    chains = np.loadtxt(chain_file)
    if chains.ndim == 1:
        chains = chains.reshape(1, -1)

    if model == 'car':
        names  = ['Omega_m', 'R_b', 'r_d', 'H0', 'S8', 'IA_bias']
        labels = [r'\Omega_m', r'R_b', r'r_d\,[{\rm Mpc}]',
                  r'H_0\,[{\rm km/s/Mpc}]', r'S_8', r'b_{\rm IA}']
        # Primary parameters are first 2; rest are derived
        samples = MCSamples(samples=chains[:, :len(names)],
                            names=names, labels=labels,
                            name_tag='CAR')
        params_to_plot = ['Omega_m', 'R_b', 'H0', 'S8']
    else:
        names  = ['Omega_b_h2', 'Omega_c_h2', 'H0', 'tau', 'A_s', 'n_s']
        labels = [r'\Omega_b h^2', r'\Omega_c h^2', r'H_0', r'\tau',
                  r'10^9 A_s', r'n_s']
        n_use  = min(len(names), chains.shape[1])
        samples = MCSamples(samples=chains[:, :n_use],
                            names=names[:n_use], labels=labels[:n_use],
                            name_tag='ΛCDM')
        params_to_plot = names[:min(4, n_use)]

    g = plots.get_subplot_plotter(width_inch=8)
    g.settings.axes_fontsize  = 11
    g.settings.legend_fontsize = 11
    g.triangle_plot([samples], params_to_plot,
                    filled=True, legend_labels=[model.upper()])

    # Add CAR prediction lines for key parameters
    if model == 'car':
        preds = CAR_predictions()
        ax_map = {name: g.get_axes_for_params(name, name)
                  for name in params_to_plot if name in params_to_plot}
        pred_vals = {
            'Omega_m': preds.get('Omega_m', 0.312),
            'R_b':     preds['R_b0'],
            'H0':      preds['H0'],
            'S8':      preds['S8'],
        }
        # Mark BBN prior on R_b
        for ax_list in g.subplots:
            for ax in (ax_list or []):
                if ax is None:
                    continue
                ax.axvline(0.260, color='red', ls=':', alpha=0.5, lw=1.2)

    import matplotlib.pyplot as plt
    plt.savefig(output_path, bbox_inches='tight', dpi=200)
    plt.close()
    print(f"  Posterior plot saved: {output_path}")


def _plot_fallback(output_path: str, model: str) -> None:
    """
    Generate a demonstration posterior plot using Gaussian approximation
    at the paper's best-fit values. Used when chain files are not available.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.stats import multivariate_normal

    preds = CAR_predictions()
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    fig.suptitle(f'{model.upper()} Posterior Distributions (demo)', fontsize=13)

    if model == 'car':
        params = [
            ('$\\Omega_m$',         0.312, 0.009,  0.28, 0.35),
            ('$R_b$',               0.260, 0.002,  0.254, 0.266),
            ('$H_0$ [km/s/Mpc]',   70.4,  0.4,    69.0, 72.0),
            ('$S_8$',               0.783, 0.015,  0.74, 0.82),
        ]
        color = '#1B6CA8'
    else:
        params = [
            ('$\\Omega_m$',        0.315, 0.007, 0.29, 0.34),
            ('$\\Omega_b h^2$',    0.0224, 0.0002, 0.021, 0.024),
            ('$H_0$ [km/s/Mpc]',  67.4,  0.5,   65.5, 69.5),
            ('$S_8$',              0.832, 0.013, 0.79, 0.87),
        ]
        color = '#888888'

    for ax, (label, mu, sigma, xlo, xhi) in zip(axes.flat, params):
        x = np.linspace(xlo, xhi, 300)
        y = np.exp(-0.5 * ((x - mu) / sigma)**2)
        ax.fill_between(x, y, alpha=0.3, color=color)
        ax.plot(x, y, color=color, lw=2)
        ax.axvline(mu, color=color, lw=1.5, ls='--')
        ax.set_xlabel(label, fontsize=12)
        ax.set_ylabel('Posterior (normalised)', fontsize=10)
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(0, 1.15)
        ax.text(0.97, 0.92, f'${mu:.3g} \\pm {sigma:.2g}$',
                transform=ax.transAxes, ha='right', fontsize=11)
        ax.grid(True, alpha=0.3)

    # Add observational constraints to H0 and S8
    for ax, (label, mu, sigma, xlo, xhi) in zip(axes.flat, params):
        if 'H_0' in label:
            ax.axvspan(73.0 - 1.0, 73.0 + 1.0, alpha=0.12, color='red',
                       label='SH0ES')
            ax.axvspan(67.4 - 0.5, 67.4 + 0.5, alpha=0.12, color='gray',
                       label='Planck ΛCDM')
            ax.legend(fontsize=9)
        if 'S_8' in label:
            ax.axvspan(0.780 - 0.012, 0.780 + 0.012, alpha=0.15, color='orange',
                       label='DES-Y6')
            ax.axvspan(0.832 - 0.013, 0.832 + 0.013, alpha=0.10, color='gray',
                       label='Planck')
            ax.legend(fontsize=9)

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=200)
    plt.close()
    print(f"  Demo posterior plot saved: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot posteriors from PolyChord chains")
    parser.add_argument("--input",  required=True,
                        help="PolyChord chain directory")
    parser.add_argument("--output", required=True,
                        help="Output PDF/PNG path")
    parser.add_argument("--model",  default="car",
                        choices=["car", "lcdm"])
    args = parser.parse_args()

    plot_car_posteriors(args.input, args.output, args.model)
