"""
data/download_all.sh → replaced by this Python script for portability.

data/mock_data_generator.py
============================
Generate mock data vectors at CAR best-fit values for pipeline testing.
Produces HDF5 files matching the format expected by each likelihood module.

Usage
-----
    python data/mock_data_generator.py --output data/
    python data/mock_data_generator.py --output data/ --seed 42 --noise 0.02

Author : DR JM NIPOK | License: GPL-3.0
"""

import argparse
import numpy as np
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sct_core import CAR_predictions

try:
    import h5py
    HDF5_AVAILABLE = True
except ImportError:
    HDF5_AVAILABLE = False
    print("Warning: h5py not installed. Saving as numpy .npz instead.")


def generate_desi_mock(params: dict, output_path: str, seed: int = 42,
                       noise_level: float = 0.02) -> None:
    """Generate mock DESI-DR2 BAO data at CAR best-fit."""
    from likelihoods.desi_dr2_bao import DESIDR2BAOLikelihood, DESI_DR2_DATA

    rng = np.random.default_rng(seed)
    lik = DESIDR2BAOLikelihood(verbose=False)
    lik.load_data()

    # Noisify the theory vector
    theory = lik.theory(params)
    noise  = noise_level * rng.standard_normal(len(theory))
    mock   = theory * (1.0 + noise)
    cov    = np.diag((noise_level * np.abs(theory))**2)

    if HDF5_AVAILABLE:
        with h5py.File(output_path, 'w') as f:
            f.create_dataset('data_vector', data=mock)
            f.create_dataset('covariance',  data=cov)
            f.attrs['n_data']      = len(mock)
            f.attrs['noise_level'] = noise_level
            f.attrs['seed']        = seed
            f.attrs['S8_input']    = params['S8']
            f.attrs['H0_input']    = params['H0']
            f.attrs['r_d_input']   = params['r_d_Mpc']
        print(f"  Written: {output_path} (n={len(mock)})")
    else:
        np.savez(output_path.replace('.hdf5','.npz'),
                 data_vector=mock, covariance=cov)
        print(f"  Written: {output_path.replace('.hdf5','.npz')} (n={len(mock)})")


def generate_des_mock(params: dict, output_path: str, seed: int = 43,
                      noise_level: float = 0.05, n_data: int = 1800) -> None:
    """Generate mock DES-Y6 3×2pt data at CAR best-fit."""
    rng = np.random.default_rng(seed)
    # Simple mock: Gaussian noise around CAR S₈-scaled signal
    S8_ratio = (params['S8'] / 0.832)**2
    base_amplitude = 1e-5 * S8_ratio
    mock = base_amplitude * np.ones(n_data) * (1.0 + noise_level * rng.standard_normal(n_data))
    cov  = np.diag((noise_level * base_amplitude * np.ones(n_data))**2)

    if HDF5_AVAILABLE:
        n_xi = n_data // 4
        with h5py.File(output_path, 'w') as f:
            f.create_dataset('xi_plus',  data=mock[:n_xi])
            f.create_dataset('xi_minus', data=mock[n_xi:2*n_xi])
            f.create_dataset('gamma_t',  data=mock[2*n_xi:3*n_xi])
            f.create_dataset('w_theta',  data=mock[3*n_xi:])
            f.create_dataset('covariance', data=cov)
            f.attrs['S8_input'] = params['S8']
            f.attrs['IA_bias']  = params['IA_bias']
        print(f"  Written: {output_path} (n={n_data})")


def print_download_instructions() -> None:
    """Print instructions for downloading real data."""
    print("""
╔══════════════════════════════════════════════════════════════════════╗
║           Real Data Download Instructions                            ║
╠══════════════════════════════════════════════════════════════════════╣
║                                                                      ║
║  DESI-DR2 BAO                                                        ║
║    URL : https://data.desi.lbl.gov/public/dr2/bao/                  ║
║    File: bao_measurements_dr2.txt                                    ║
║    DOI : arXiv:2603.12458                                            ║
║                                                                      ║
║  DES-Y6 3×2pt                                                        ║
║    URL : https://des.ncsa.illinois.edu/releases/y6                   ║
║    File: des_y6_3x2pt_data_vector.fits                               ║
║    DOI : arXiv:2603.08941                                            ║
║                                                                      ║
║  HSC-Y3 Weak Lensing                                                 ║
║    URL : https://hsc-release.mtk.nao.ac.jp/dataservice/y3/wl/       ║
║    File: hsc_y3_cosmic_shear.hdf5                                    ║
║    DOI : ApJ 960, 62 (2025)                                          ║
║                                                                      ║
║  KiDS-DR5 Weak Lensing                                               ║
║    URL : https://kids.strw.leidenuniv.nl/DR5/                        ║
║    File: KiDS-1000_Cosmology_xi_pm.fits                              ║
║    DOI : A&A 688, A120 (2026)                                        ║
║                                                                      ║
║  Planck PR4 (NPIPE)                                                  ║
║    URL : https://pla.esac.esa.int (requires ESA login)               ║
║    File: COM_Likelihood_Data-PR4_2020.tar.gz                         ║
║    DOI : A&A 687, A160 (2024)                                        ║
║                                                                      ║
║  After downloading, convert to HDF5 format using:                   ║
║    python data/convert_to_hdf5.py --survey desi --input <file>      ║
╚══════════════════════════════════════════════════════════════════════╝
""")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate mock data for testing")
    parser.add_argument("--output",      default="data/",
                        help="Output directory for mock HDF5 files")
    parser.add_argument("--seed",        type=int,   default=42)
    parser.add_argument("--noise",       type=float, default=0.02,
                        help="Relative noise level (default: 0.02 = 2%%)")
    parser.add_argument("--real-data",   action="store_true",
                        help="Print instructions for downloading real data")
    args = parser.parse_args()

    if args.real_data:
        print_download_instructions()
        sys.exit(0)

    os.makedirs(args.output, exist_ok=True)
    params = CAR_predictions()

    print(f"Generating mock data at CAR best-fit:")
    print(f"  H₀ = {params['H0']:.2f}, S₈ = {params['S8']:.3f}, r_d = {params['r_d_Mpc']:.2f} Mpc")
    print()

    generate_desi_mock(params,
                       os.path.join(args.output, 'desi_dr2_bao.hdf5'),
                       seed=args.seed, noise_level=args.noise)
    generate_des_mock(params,
                      os.path.join(args.output, 'des_y6_3x2pt.hdf5'),
                      seed=args.seed+1, noise_level=args.noise*2)

    print()
    print("Mock data generation complete.")
    print("For real data, run: python data/mock_data_generator.py --real-data")
