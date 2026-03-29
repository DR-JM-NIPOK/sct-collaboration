"""
figures/make_all_figures.py
============================
Reproduce all paper figures for SCT Cosmology Series Paper #16.

Usage
-----
    python figures/make_all_figures.py            # all figures
    python figures/make_all_figures.py --fig 1    # specific figure
    python figures/make_all_figures.py --format pdf --dpi 300

Figures produced
----------------
    fig1_tension_summary.pdf   — H₀, S₈, r_d tension triangle
    fig2_sound_horizon.pdf     — c_s(z) and r_d integrand comparison
    fig3_posteriors.pdf        — CAR vs ΛCDM posterior contours
    fig4_residuals.pdf         — Standardised residuals by dataset
    fig5_forecasts.pdf         — Euclid/LSST/CMB-S4 kill switch forecasts

Author : DR JM NIPOK | License: GPL-3.0
"""

import argparse
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from sct_core import CAR_predictions, lcdm_reference, cs2_CAR as cs_squared, R_b_of_z

# ─── Style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':    'serif',
    'font.size':      12,
    'axes.labelsize': 13,
    'legend.fontsize':11,
    'figure.dpi':     150,
    'lines.linewidth':1.8,
})

CAR_COLOR  = '#1B6CA8'   # blue
LCDM_COLOR = '#888888'   # grey
OBS_COLOR  = '#CC3311'   # red


def fig1_tension_summary(output_dir: str, fmt: str = 'pdf') -> None:
    """
    Figure 1: H₀, S₈, r_d tension summary.
    Shows CAR predictions vs SH0ES, Planck ΛCDM, and survey data.
    """
    preds = CAR_predictions()
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    fig.suptitle('CAR Parameter-Free Predictions vs Observations', fontsize=13, y=1.02)

    # Panel 1: H₀
    ax = axes[0]
    data_pts = [
        ('SH0ES',      73.0, 1.0),
        ('Planck ΛCDM',67.4, 0.5),
        ('DESI-DR2',   68.5, 1.2),
    ]
    ys = np.arange(len(data_pts))
    for i, (label, val, err) in enumerate(data_pts):
        ax.errorbar(val, i, xerr=err, fmt='o', color=OBS_COLOR, ms=7, capsize=4)
        ax.text(val + err + 0.1, i, label, va='center', fontsize=10)
    ax.axvline(preds['H0'], color=CAR_COLOR, lw=2.5, ls='--', label=f"CAR: {preds['H0']:.1f}")
    ax.set_xlabel(r'$H_0$ [km/s/Mpc]')
    ax.set_yticks([]); ax.legend(fontsize=10)
    ax.set_xlim(64, 76)
    ax.set_title(r'$H_0$ tension', fontsize=12)

    # Panel 2: S₈
    ax = axes[1]
    data_pts_s8 = [
        ('Planck ΛCDM', 0.832, 0.013),
        ('DES-Y6',      0.780, 0.012),
        ('HSC-Y3',      0.776, 0.032),
        ('KiDS-DR5',    0.788, 0.014),
    ]
    for i, (label, val, err) in enumerate(data_pts_s8):
        color = OBS_COLOR if 'ΛCDM' not in label else LCDM_COLOR
        ax.errorbar(val, i, xerr=err, fmt='s', color=color, ms=7, capsize=4)
        ax.text(val + err + 0.002, i, label, va='center', fontsize=10)
    ax.axvline(preds['S8'], color=CAR_COLOR, lw=2.5, ls='--', label=f"CAR: {preds['S8']:.3f}")
    ax.set_xlabel(r'$S_8 \equiv \sigma_8\sqrt{\Omega_m/0.3}$')
    ax.set_yticks([]); ax.legend(fontsize=10)
    ax.set_xlim(0.74, 0.87)
    ax.set_title(r'$S_8$ tension', fontsize=12)

    # Panel 3: r_d
    ax = axes[2]
    data_pts_rd = [
        ('Planck ΛCDM', 150.0, 0.4),
        ('DESI-DR2',    147.0, 1.0),
    ]
    for i, (label, val, err) in enumerate(data_pts_rd):
        color = LCDM_COLOR if 'ΛCDM' in label else OBS_COLOR
        ax.errorbar(val, i, xerr=err, fmt='^', color=color, ms=8, capsize=4)
        ax.text(val + err + 0.05, i, label, va='center', fontsize=10)
    ax.axvline(preds['r_d_Mpc'], color=CAR_COLOR, lw=2.5, ls='--',
               label=f"CAR: {preds['r_d_Mpc']:.1f}")
    ax.set_xlabel(r'$r_d$ [Mpc]')
    ax.set_yticks([]); ax.legend(fontsize=10)
    ax.set_xlim(144, 152)
    ax.set_title(r'$r_d$ tension', fontsize=12)

    plt.tight_layout()
    out = os.path.join(output_dir, f'fig1_tension_summary.{fmt}')
    plt.savefig(out, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out}")


def fig2_sound_horizon(output_dir: str, fmt: str = 'pdf') -> None:
    """
    Figure 2: Sound speed c_s(z) and integrand c_s(z)/H(z).
    Compares CAR vs standard ΛCDM.
    """
    from sct_core import E_of_z as hubble_factor, BBN_OMEGA_B_H2, PLANCK_OMEGA_GAM_H2, PLANCK_OMEGA_M

    R_b0 = (4 * BBN_OMEGA_B_H2) / (3 * PLANCK_OMEGA_GAM_H2)
    h = 0.704
    Omega_r = (PLANCK_OMEGA_GAM_H2 / h**2) * (1 + 0.2271 * 3.044)
    z_arr = np.logspace(-0.5, 3.1, 500)

    cs2_car  = np.array([cs_squared(R_b0, z) for z in z_arr])
    cs2_lcdm = np.array([1.0 / (3.0 * (1.0 + R_b_of_z(R_b0, z))) for z in z_arr])
    Hz       = np.array([hubble_factor(z, PLANCK_OMEGA_M, Omega_r) for z in z_arr])

    integrand_car  = np.sqrt(cs2_car)  / Hz
    integrand_lcdm = np.sqrt(cs2_lcdm) / Hz

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    ax1.semilogx(z_arr, np.sqrt(cs2_car),  color=CAR_COLOR,  label='CAR: $c_s = \\sqrt{(1+R_b)/3}$')
    ax1.semilogx(z_arr, np.sqrt(cs2_lcdm), color=LCDM_COLOR, ls='--', label='ΛCDM: $c_s = 1/\\sqrt{3(1+R_b)}$')
    ax1.axvline(1089, color='gray', ls=':', alpha=0.7, label='$z_*=1089$')
    ax1.set_xlabel('Redshift $z$')
    ax1.set_ylabel('Sound speed $c_s$ $[c]$')
    ax1.set_title('CAR vs ΛCDM sound speed')
    ax1.legend(); ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0.3, 1300)
    ax1.set_ylim(0.0, 0.75)

    ax2.semilogx(z_arr, integrand_car,  color=CAR_COLOR,  label='CAR integrand')
    ax2.semilogx(z_arr, integrand_lcdm, color=LCDM_COLOR, ls='--', label='ΛCDM integrand')
    # Mark z* and show cumulative contribution
    idx_zstar = np.argmin(np.abs(z_arr - 1089))
    from scipy.integrate import cumulative_trapezoid
    cum_car = cumulative_trapezoid(integrand_car[::-1], z_arr[::-1], initial=0)[::-1]
    ax2.fill_between(z_arr, integrand_car, alpha=0.15, color=CAR_COLOR)
    ax2.axvline(1089, color='gray', ls=':', alpha=0.7)
    ax2.set_xlabel('Redshift $z$')
    ax2.set_ylabel('$c_s(z)/H(z)$ $[c/H_0]$')
    ax2.set_title(f'Sound horizon integrand\n(shaded area = $r_d$)')
    ax2.legend(); ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0.3, 1300)

    plt.tight_layout()
    out = os.path.join(output_dir, f'fig2_sound_horizon.{fmt}')
    plt.savefig(out, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out}")


def fig5_forecasts(output_dir: str, fmt: str = 'pdf') -> None:
    """
    Figure 5: Falsifiability forecasts for Euclid, LSST, CMB-S4.
    """
    preds = CAR_predictions()
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    fig.suptitle('CAR Falsifiable Forecasts', fontsize=13)

    # Panel 1: Euclid Y1 S₈
    ax = axes[0]
    S8_pred   = 0.781
    S8_sigma  = 0.009
    S8_kill   = 0.81
    x = np.linspace(0.74, 0.86, 300)
    y = np.exp(-0.5*((x - S8_pred)/S8_sigma)**2)
    ax.fill_between(x, y, alpha=0.3, color=CAR_COLOR, label=f'CAR: {S8_pred} ± {S8_sigma}')
    ax.plot(x, y, color=CAR_COLOR, lw=2)
    ax.axvline(S8_kill, color='red', ls='--', lw=2, label=f'Kill switch: $S_8 > {S8_kill}$')
    ax.fill_betweenx([0, 1.1], S8_kill, 0.86, alpha=0.15, color='red')
    ax.set_xlabel(r'$S_8$'); ax.set_ylabel('Posterior (normalised)')
    ax.set_title('Euclid Year 1 (2027)'); ax.legend(fontsize=10)
    ax.set_ylim(0, 1.15); ax.set_xlim(0.74, 0.86)

    # Panel 2: LSST Y1
    ax = axes[1]
    S8_lsst  = 0.780
    sig_lsst = 0.006
    x2 = np.linspace(0.75, 0.82, 300)
    y2 = np.exp(-0.5*((x2 - S8_lsst)/sig_lsst)**2)
    ax.fill_between(x2, y2, alpha=0.3, color=CAR_COLOR, label=f'CAR: {S8_lsst} ± {sig_lsst}')
    ax.plot(x2, y2, color=CAR_COLOR, lw=2)
    ax.axvline(0.805, color='red', ls='--', lw=2, label='$\\chi^2/dof > 1.05$')
    ax.set_xlabel(r'$S_8$')
    ax.set_title('LSST Year 1 (2028)'); ax.legend(fontsize=10)
    ax.set_ylim(0, 1.15)

    # Panel 3: CMB-S4 r_d
    ax = axes[2]
    rd_pred  = 149.2
    rd_sigma = 0.2
    rd_kill_lo, rd_kill_hi = 148.5, 150.0
    x3 = np.linspace(147.5, 151.0, 300)
    y3 = np.exp(-0.5*((x3 - rd_pred)/rd_sigma)**2)
    ax.fill_between(x3, y3, alpha=0.3, color=CAR_COLOR, label=f'CAR: {rd_pred} ± {rd_sigma}')
    ax.plot(x3, y3, color=CAR_COLOR, lw=2)
    ax.axvline(rd_kill_lo, color='red', ls='--', lw=2, label='Kill region')
    ax.axvline(rd_kill_hi, color='red', ls='--', lw=2)
    ax.fill_betweenx([0,1.1], 147.5, rd_kill_lo, alpha=0.15, color='red')
    ax.fill_betweenx([0,1.1], rd_kill_hi, 151.0, alpha=0.15, color='red')
    ax.set_xlabel(r'$r_d$ [Mpc]')
    ax.set_title('CMB-S4 (2030)'); ax.legend(fontsize=10)
    ax.set_ylim(0, 1.15)

    plt.tight_layout()
    out = os.path.join(output_dir, f'fig5_forecasts.{fmt}')
    plt.savefig(out, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out}")


FIGURE_FUNCTIONS = {
    1: fig1_tension_summary,
    2: fig2_sound_horizon,
    5: fig5_forecasts,
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate all paper figures")
    parser.add_argument("--output", default="output/figures",
                        help="Output directory for figures")
    parser.add_argument("--format", default="pdf", choices=["pdf", "png", "svg"],
                        help="Output format")
    parser.add_argument("--fig",    type=int, default=None,
                        help="Generate specific figure only (1-5)")
    parser.add_argument("--dpi",    type=int, default=150)
    args = parser.parse_args()

    plt.rcParams['figure.dpi'] = args.dpi
    os.makedirs(args.output, exist_ok=True)

    figs_to_run = [args.fig] if args.fig else sorted(FIGURE_FUNCTIONS.keys())

    for fig_num in figs_to_run:
        if fig_num in FIGURE_FUNCTIONS:
            print(f"Generating Figure {fig_num}...")
            FIGURE_FUNCTIONS[fig_num](args.output, args.format)
        else:
            print(f"Figure {fig_num} not yet implemented (requires full chain output).")

    print(f"\nDone. Figures saved to {args.output}/")
