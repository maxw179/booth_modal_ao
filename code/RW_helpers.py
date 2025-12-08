from multiprocessing import Pool, cpu_count, get_context
import numpy as np
from scipy.special import j0, j1, jn


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
    r_bfp =  f * n * np.sin(theta)  # mm
    #rr_bfp = f * np.sin(theta)  # mm
    return np.exp(-4 * np.log(2) * (r_bfp**2) / fwhm_pupil**2)

#r_p cos(epsilon)
def aberration_free_phase(k, x, y, z, theta, phi):
    return \
        k * (x * np.sin(theta) * np.cos(phi)
        + y * np.sin(theta) * np.sin(phi)
        + z * np.cos(theta))

def E_integrate(x,y,z, alpha, k, f, n, fwhm_pupil,
                n_theta = 128, n_phi = 128, aberration_map = None):
    #define the grid of allowable thetas and phis
    theta = np.linspace(0, alpha, n_theta)
    phi = np.linspace(0, 2 * np.pi, n_phi)

    #get the riemann integration element
    d_theta = theta[1] - theta[0]
    d_phi = phi[1] - phi[0]

    #get the pairs of theta, phi to integrate over
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing = 'ij')
    #helpful trig things for later
    cosT = np.cos(theta_grid)
    sinT = np.sin(theta_grid)
    #get the factors out in front
    pre_factor = -1j * k * f / (2 * np.pi)
    inside_factor = gaussian_amplitude(theta_grid, f, n, fwhm_pupil) * np.sqrt(cosT) * sinT
    #check if there is extra aberration to add
    aberration_phase = 0.0
    if not(aberration_map is None):
        aberration_phase = aberration_map(theta_grid, phi_grid)
    #calculate the total phase offset:
    default_phase = aberration_free_phase(k, x, y, z, theta_grid, phi_grid)
    phase = np.exp(1j * (default_phase + aberration_phase))

    #get the strength factors
    a_x, a_y, a_z = strength_angular(theta_grid, phi_grid)
    #take the integral!
    E_x = pre_factor * np.sum(inside_factor * a_x * phase) * d_theta * d_phi
    E_y = pre_factor * np.sum(inside_factor * a_y * phase) * d_theta * d_phi
    E_z = pre_factor * np.sum(inside_factor * a_z * phase) * d_theta * d_phi
    return E_x, E_y, E_z


def intensity_grid(L_ffp, grid_ffp, alpha, k, f, n, fwhm_pupil, theta_grid_size, N_order, aberration_map = None):
    x = np.linspace(-L_ffp/2, L_ffp/2, grid_ffp)
    y = np.linspace(-L_ffp/2, L_ffp/2, grid_ffp)
    intensity_map = np.zeros((len(x), len(y)))

    for i in range(len(x)):
        for j in range(len(y)):
            x_p = x[i]
            y_p = y[j]
            #get the integrated electric field
            E_x, E_y, E_z = E_integrate(x_p,y_p,0,alpha,k,f,n,fwhm_pupil, theta_grid_size, theta_grid_size, aberration_map)
            I1 = np.abs(E_x)**2 + np.abs(E_y)**2 + np.abs(E_z)**2
            I_np = I1**N_order
            intensity_map[i, j] = I_np

    return x, y, intensity_map.T

def _init_worker(aberration_map):
    global _aberration_map
    _aberration_map = aberration_map


def row_intensity_helper(args):
    (x_i, y_array, z,
     alpha, k, f, n, fwhm_pupil,
     theta_grid_size, N_order) = args

    row = np.empty_like(y_array, dtype=float)

    for idx, y_j in enumerate(y_array):
        E_x, E_y, E_z = E_integrate(
            x_i, y_j, z,
            alpha, k, f, n,
            fwhm_pupil,
            theta_grid_size, theta_grid_size,
            _aberration_map
        )
        I1 = np.abs(E_x)**2 + np.abs(E_y)**2 + np.abs(E_z)**2
        row[idx] = I1**N_order

    return row

def intensity_grid_parallel(
    L_ffp, grid_ffp,
    alpha, k, f, n, fwhm_pupil,
    theta_grid_size, N_order,
    aberration_map=None,
    n_procs=None,
):
    x = np.linspace(-L_ffp / 2, L_ffp / 2, grid_ffp)
    y = np.linspace(-L_ffp / 2, L_ffp / 2, grid_ffp)
    z = 0.0

    # One task per x-row
    tasks = [
        (x_i, y, z,
            alpha, k, f, n, fwhm_pupil,
            theta_grid_size, N_order)
        for x_i in x
    ]

    if n_procs is None:
        n_procs = cpu_count()

    ctx = get_context("spawn")
    with ctx.Pool(
        processes=n_procs,
        initializer=_init_worker,
        initargs=(aberration_map,)
    ) as pool:
        intensities_rows = pool.map(row_intensity_helper, tasks)

    #stack rows into a 2D array: shape (len(x), len(y))
    intensity_map = np.vstack(intensities_rows)

    return x, y, intensity_map.T