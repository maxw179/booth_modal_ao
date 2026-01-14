from multiprocessing import Pool, cpu_count, get_context
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j0, j1, jn
import subprocess
import math
import Zernike_helpers

def strength_angular(theta, phi):
    cosT = np.cos(theta)
    sinT = np.sin(theta)
    cosP = np.cos(phi)
    sinP = np.sin(phi)

    a_x = cosT + (1 - cosT) * sinP**2
    a_y = (cosT - 1) * cosP * sinP
    a_z = -sinT * cosP

    return a_x, a_y, a_z

def gaussian_amplitude(theta, f, n, fwhm_pupil):
    r_bfp = f * n * np.sin(theta)  # mm
    return np.exp(-4 * np.log(2) * (r_bfp**2) / fwhm_pupil**2)

#r_p cos(epsilon)
def aberration_free_phase(k, x, y, z, theta, phi):
    return \
        k * (x * np.sin(theta) * np.cos(phi)
        + y * np.sin(theta) * np.sin(phi)
        + z * np.cos(theta))

def make_pupil_cache(alpha, k, f, n, fwhm_pupil,
                     n_theta=128, n_phi=128,
                     aberration_map=None):
    #define the grid of allowable thetas and phis
    theta = np.linspace(0, alpha, n_theta, endpoint=False)
    phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)

    #get the riemann integration element
    d_theta = theta[1] - theta[0]
    d_phi = phi[1] - phi[0]

    #get the pairs of theta, phi to integrate over
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")

    #helpful trig things for later
    cosT = np.cos(theta_grid)
    sinT = np.sin(theta_grid)

    #get the factors out in front
    pre_factor = -1j * k * f / (2 * np.pi)
    inside_factor = gaussian_amplitude(theta_grid, f, n, fwhm_pupil) * np.sqrt(cosT) * sinT

    #precompute strength factors
    a_x, a_y, a_z = strength_angular(theta_grid, phi_grid)

    #precompute aberration phase on the pupil grid
    aberration_phase = 0.0
    if aberration_map is not None:
        aberration_phase = aberration_map(theta_grid, phi_grid)


    return {
        "theta_grid": theta_grid,
        "phi_grid": phi_grid,
        "d_theta": d_theta,
        "d_phi": d_phi,
        "pre_factor": pre_factor,
        "inside_factor": inside_factor,
        "a_x": a_x,
        "a_y": a_y,
        "a_z": a_z,
        "aberration_phase": aberration_phase,
    }


def E_integrate_cached(x, y, z, k, cache):
    
    theta_grid = cache["theta_grid"]
    phi_grid = cache["phi_grid"]

    default_phase = aberration_free_phase(k, x, y, z, theta_grid, phi_grid)
    phase = np.exp(1j * (default_phase + cache["aberration_phase"]))

    # take the integral
    E_x = cache["pre_factor"] * np.sum(cache["inside_factor"] * cache["a_x"] * phase) * cache["d_theta"] * cache["d_phi"]
    E_y = cache["pre_factor"] * np.sum(cache["inside_factor"] * cache["a_y"] * phase) * cache["d_theta"] * cache["d_phi"]
    E_z = cache["pre_factor"] * np.sum(cache["inside_factor"] * cache["a_z"] * phase) * cache["d_theta"] * cache["d_phi"]
    return E_x, E_y, E_z


def intensity_grid(L_ffp, grid_ffp, alpha, k, f, n, fwhm_pupil, theta_grid_size, N_order, aberration_map=None):
    x = np.linspace(-L_ffp / 2, L_ffp / 2, grid_ffp)
    y = np.linspace(-L_ffp / 2, L_ffp / 2, grid_ffp)
    intensity_map = np.zeros((len(x), len(y)))

    #build cache once for the worker
    cache = make_pupil_cache(
        alpha=alpha, k=k, f=f, n=n, fwhm_pupil=fwhm_pupil,
        n_theta=theta_grid_size, n_phi=theta_grid_size,
        aberration_map=aberration_map
    )

    for i in range(len(x)):
        for j in range(len(y)):
            x_p = x[i]
            y_p = y[j]
            #get the integrated electric field (with cached pupil)
            E_x, E_y, E_z = E_integrate_cached(x_p, y_p, 0.0, k, cache)
            I1 = np.abs(E_x) ** 2 + np.abs(E_y) ** 2 + np.abs(E_z) ** 2
            I_np = I1 ** N_order
            intensity_map[i, j] = I_np

    return x, y, np.flip(intensity_map.T)

def _init_worker(aberration_kind, params, alpha, k, f, n, fwhm_pupil, theta_grid_size):
    global _aberration_map, _pupil_cache

    if aberration_kind == "Zernike":
        _aberration_map = Zernike_helpers.create_zernike_function(*params)
    elif aberration_kind == "Random":
        _aberration_map = Zernike_helpers.generate_random_aberration(*params)
    else:
        _aberration_map = None

    #precompute cache for each worker
    _pupil_cache = make_pupil_cache(
        alpha=alpha, k=k, f=f, n=n, fwhm_pupil=fwhm_pupil,
        n_theta=theta_grid_size, n_phi=theta_grid_size,
        aberration_map=_aberration_map
    )

def row_intensity_helper(args):
    (x_i, y_array, z, k, N_order) = args

    row = np.empty_like(y_array, dtype=float)

    for idx, y_j in enumerate(y_array):
        E_x, E_y, E_z = E_integrate_cached(x_i, y_j, z, k, _pupil_cache)
        I1 = np.abs(E_x) ** 2 + np.abs(E_y) ** 2 + np.abs(E_z) ** 2
        row[idx] = I1 ** N_order

    return row

def intensity_grid_parallel(
    L_ffp, grid_ffp,
    alpha, k, f, n, fwhm_pupil,
    theta_grid_size, N_order,
    aberration_kind=None,
    n_procs=None,
    params=None):

    #we will always need this to construct a zernike map
    if params is not None:
        params.append(alpha)

    x = np.linspace(-L_ffp / 2, L_ffp / 2, grid_ffp)
    y = np.linspace(-L_ffp / 2, L_ffp / 2, grid_ffp)
    z = 0.0

    #one task per row
    tasks = [
        (x_i, y, z, k, N_order)
        for x_i in x
    ]

    if n_procs is None:
        n_procs = cpu_count()

    ctx = get_context("spawn")
    with ctx.Pool(
        processes=n_procs,
        initializer=_init_worker,
        initargs=(aberration_kind, params, alpha, k, f, n, fwhm_pupil, theta_grid_size)
    ) as pool:
        intensities_rows = pool.map(row_intensity_helper, tasks)

    # tack rows into a 2D array: shape (len(x), len(y))
    intensity_map = np.vstack(intensities_rows)

    return x, y, np.flip(intensity_map.T)

def parallel_grid_wrapper(L_ffp, grid_ffp, alpha, k, f, n, fwhm_pupil, theta_grid_size, N_order,
                          aberration_kind=None, output="psf_output.npz", params=None, python_executable="python",
                          script_path="RW_run_parallel.py"):
    cmd = [
        python_executable, script_path,
        "--L-ffp", str(L_ffp),
        "--grid-ffp", str(grid_ffp),
        "--alpha", str(alpha),
        "--k", str(k),
        "--f", str(f),
        "--n", str(n),
        "--fwhm-pupil", str(fwhm_pupil),
        "--theta-grid-size", str(theta_grid_size),
        "--N-order", str(N_order),
        "--output", output,
        "--params", str(params),
    ]

    if aberration_kind is not None:
        cmd += ["--aberration-kind", str(aberration_kind)]

    subprocess.run(cmd, check=True)

    output_files = np.load(output)
    return output_files["x"], output_files["y"], output_files["I"]