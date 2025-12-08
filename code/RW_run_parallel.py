#!/usr/bin/env python
import argparse
import numpy as np
from RW_helpers import intensity_grid_parallel


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute PSF intensity grid using multiprocessing."
    )

    #basic grid / region params
    parser.add_argument(
        "--L-ffp",
        type=float,
        required=True,
        help="Field of view in the Fourier/focal plane (same units as x,y).",
    )
    parser.add_argument(
        "--grid-ffp",
        type=int,
        required=True,
        help="Number of points along each axis in the x-y grid.",
    )

    #optical parameters
    parser.add_argument(
        "--alpha",
        type=float,
        required=True,
        help="Half-angle of the objective (in radians).",
    )
    parser.add_argument(
        "--k",
        type=float,
        required=True,
        help="Wavenumber k = 2πn/λ (in 1/length).",
    )
    parser.add_argument(
        "--f",
        type=float,
        required=True,
        help="Focal length of the objective (same length units as x,y).",
    )
    parser.add_argument(
        "--n",
        type=float,
        required=True,
        help="Refractive index of the medium.",
    )
    parser.add_argument(
        "--fwhm-pupil",
        type=float,
        required=True,
        help="FWHM of the Gaussian at the pupil (in pupil units).",
    )

    # Integration grid + nonlinearity order
    parser.add_argument(
        "--theta-grid-size",
        type=int,
        required=True,
        help="Number of theta (and phi) samples for the angular integration.",
    )
    parser.add_argument(
        "--N-order",
        type=int,
        required=True,
        help="Order of nonlinearity (2 for 2P, 3 for 3P, etc.).",
    )

    # Optional aberration map
    parser.add_argument(
        "--aberration-map",
        type=str,
        default=None,
        help="Path to .npy file containing aberration map (or omit for none).",
    )

    # Misc
    parser.add_argument(
        "--n-procs",
        type=int,
        default=None,
        help="Number of processes to use (default: cpu_count).",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="psf_result.npz",
        help="Output .npz file to write x, y, I to (default: psf_result.npz).",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    #load aberration map if provided
    if args.aberration_map is not None:
        aberration_map = np.load(args.aberration_map)
    else:
        aberration_map = None

    #call the parallel PSF computation
    x, y, I = intensity_grid_parallel(
        L_ffp=args.L_ffp,
        grid_ffp=args.grid_ffp,
        alpha=args.alpha,
        k=args.k,
        f=args.f,
        n=args.n,
        fwhm_pupil=args.fwhm_pupil,
        theta_grid_size=args.theta_grid_size,
        N_order=args.N_order,
        aberration_map=aberration_map,
        n_procs=args.n_procs,
    )

    #save to disk for the notebook (or anything) to load later
    np.savez(args.output, x=x, y=y, I=I)

if __name__ == "__main__":
    main()
