import numpy as np 
import scipy.optimize as opt
from Zernike_helpers import *
from RW_helpers import * 
from Img_helpers import *
from Plot_helpers import *

class Aberration():
    def __init__(self, modes, strengths, alpha):
        self.modes = modes 
        self.strengths = strengths
        self.alpha = alpha
    def __str__(self):
        s = ""
        for i, mode in enumerate(self.modes):
            s += f"m={mode[0]}, n={mode[1]}: {np.round(self.strengths[i], 3)}"
        return s
    def __add__(self, other):
        return Aberration(self.modes, self.strengths + other.strengths, self.alpha)

    def construct_map(self):
        return create_zernike_function(self.modes, self.strengths, self.alpha)
    
class Microscope():
    def __init__(self, N_order, lambd, n, num_apt, length, mag, fwhm_pupil,
                 L_ffp, res, theta_grid_size):
        self.N_order = N_order 
        self.lambd =lambd
        self.n = n
        self.num_apt = num_apt 
        self.length = length 
        self.mag = mag
        self.fwhm_pupil = fwhm_pupil

        self.k = (2*n*np.pi)/lambd 
        self.f = length/mag 
        self.na_eff = num_apt 
        self.alpha = np.arcsin(self.na_eff / n)
        
        self.L_ffp = L_ffp 
        self.res = res
        self.grid_ffp = int(L_ffp/res + 1)
        self.theta_grid_size = theta_grid_size

class SLM():
    def __init__(self, microscope, raw_image_x, raw_image_y, raw_image, modes, strengths):
        self.microscope = microscope 
        self.innate_aberration = Aberration(modes, strengths, self.microscope.alpha)
        self.corrected_aberration = Aberration(modes, np.zeros_like(strengths), self.microscope.alpha)

        self.raw_image_x = raw_image_x
        self.raw_image_y = raw_image_y
        self.raw_image = raw_image

        self.psf_x, self.psf_y, self.psf_0 = self.get_uncorrected_psf(return_dims = True)
        self.img_0 = self.get_uncorrected_img()
        self.quality_0 = quality(self.img_0)
        self.ideal_quality = self._get_ideal_quality()

    def reset_correction(self):
        self.corrected_aberration.strengths = np.zeros_like(self.corrected_aberration.strengths)

    def _compute_psf(self, modes, strengths, return_dims = False):
        x, y, psf = parallel_grid_wrapper(self.microscope.L_ffp, 
                                              self.microscope.grid_ffp, 
                                              self.microscope.alpha, 
                                              self.microscope.k, 
                                              self.microscope.f, 
                                              self.microscope.n, 
                                              self.microscope.fwhm_pupil, 
                                              self.microscope.theta_grid_size, 
                                              self.microscope.N_order, "Zernike", "parallel_output/psf_output.npz",
                                              [modes.tolist(), strengths.tolist()])
        if return_dims:
            return x, y, psf
        else:
            return psf
        
    def get_uncorrected_psf(self, return_dims = False):
        return self._compute_psf(modes = self.innate_aberration.modes,
                                 strengths = self.innate_aberration.strengths,
                                 return_dims = return_dims)

    def get_uncorrected_img(self):
        return psf_convolve(self.raw_image, self.psf_0)

    def get_uncorrected_quality(self):
        return quality(self.img_0)
    
    def get_corrected_psf(self):
        return self._compute_psf(modes = self.innate_aberration.modes,
                                 strengths = self.innate_aberration.strengths + self.corrected_aberration.strengths)
    
    def get_corrected_img(self):
        corrected_psf = self.get_corrected_psf()
        return psf_convolve(self.raw_image, corrected_psf)
    
    def get_corrected_quality(self):
        return quality(self.get_corrected_img())
    
    def get_biased_psf(self, biases):
        return self._compute_psf(modes = self.innate_aberration.modes,
                                 strengths = self.innate_aberration.strengths + biases)

    def get_biased_img(self, biases):
        biased_psf = self.get_biased_psf(biases)
        return psf_convolve(self.raw_image, biased_psf)
    
    def get_biased_quality(self, biases):
        return quality(self.get_biased_img(biases))
    
    def _get_ideal_quality(self):
        perfect_psf = self._compute_psf(modes = np.array([]), strengths = np.array([]))
        perfect_img = psf_convolve(self.raw_image, perfect_psf)
        return quality(perfect_img)

    def run_booth(self, num_iterations, corrected_modes = None, init_bias = 0.2, zeta = 0.25):
        print("===== Running Booth 2002 Modal AO Algorithm =====")
        if corrected_modes is None:
            print("No specific modes inputted. Assuming all modes should be corrected.")
            corrected_modes = self.innate_aberration.modes
        else:
            print("Noted that only some modes should be corrected.")
        
        print("Original Aberration:")
        print(self.innate_aberration)
        #each subarray is one iteration worth of adjustments to be made
        new_corrections = np.zeros((num_iterations, len(self.corrected_aberration.strengths)))
        
        #data on the slm signal applied and their corresponding metrics 
            #each iteration, then all the modes, then the measurement
        slm_output_data = np.zeros((num_iterations + 1, len(self.innate_aberration.modes), 3))
        quality_data = np.zeros((num_iterations + 1, len(self.innate_aberration.modes), 3))
        
        bias = init_bias
        for iteration in range(num_iterations):
            print(f"===== Iteration {iteration + 1}=====")
            #get unbiased image
            print("Getting unbiased image at: " + str(self.corrected_aberration))
            M_0 = self.get_corrected_quality()
            #put the data in the arrays
            slm_output_data[iteration, :, 0] = self.corrected_aberration.strengths 
            quality_data[iteration, :, 0] = M_0
            #print the image quality for the user
            print("Image Quality:", np.round(M_0/self.quality_0, 3))
            #optimize each mode separately
            print(f"Applying biases: {np.round(bias, 3)} waves")
            for mode_index, mode in enumerate(self.innate_aberration.modes):
                if mode.tolist() in corrected_modes:
                    print(f"\t===Optimizing m={mode[0]}, n={mode[1]}===")
                    #initialize the arrays with the biases
                    plus_bias = np.copy(self.corrected_aberration.strengths)
                    minus_bias = np.copy(self.corrected_aberration.strengths)
                    plus_bias[mode_index] += bias 
                    minus_bias[mode_index] -= bias
                    print("\tExploring plus bias...")
                    M_plus =  self.get_biased_quality(plus_bias)
                    slm_output_data[iteration, mode_index, 1] = plus_bias[mode_index]
                    quality_data[iteration, mode_index, 1] = M_plus
                    print("\tImage Quality:", np.round(M_plus/self.quality_0, 3))
                    print("\tExploring minus bias...")
                    M_minus =  self.get_biased_quality(minus_bias)
                    slm_output_data[iteration, mode_index, 2] = minus_bias[mode_index]
                    quality_data[iteration, mode_index, 2] = M_minus
                    print("\tImage Quality:", np.round(M_minus/self.quality_0, 3))
                    print("\tEvaluating next step...")
                    step = (M_plus - M_minus)/(self.quality_0*2*bias) * zeta
                    #don't step by more than one half
                    if step > 0.5:
                        step = 0.5
                    if step < -0.5:
                        step = -0.5
                    print(f"\tStepping by {np.round(step, 3)} to {np.round(self.corrected_aberration.strengths[mode_index] + step, 3)}")
                    new_corrections[iteration, mode_index] = self.corrected_aberration.strengths[mode_index] + step
                    print("Applying corrections...")
                    self.corrected_aberration.strengths[mode_index] = new_corrections[iteration, mode_index]
            bias = np.max([0.05, bias * (4/5)])
            
        print("Getting final image at: " + str(self.corrected_aberration))
        M_f = self.get_corrected_quality()
        slm_output_data[-1, :, 0] = self.corrected_aberration.strengths 
        quality_data[-1, :, 0] = M_f
        print("Final image quality:", np.round(M_f/self.quality_0, 3))
        return slm_output_data, quality_data, new_corrections
    
    def get_corrected_map(self):
        return self.corrected_aberration.construct_map()
    
    def get_uncorrected_map(self):
        return self.innate_aberration.construct_map()
    
    def get_empty_map(self):
        return create_zernike_function([], [], self.microscope.alpha)
    
    def zernike_RMS(z_1, z_2, alpha, n=100):
        x = np.linspace(-1, 1, n)
        y = np.linspace(-1, 1, n)

        x_grid, y_grid = np.meshgrid(x, y, indexing="xy")
        rho_grid = np.sqrt(x_grid**2 + y_grid**2)
        phi_grid = np.arctan2(y_grid, x_grid)

        mask = rho_grid <= 1.0
        theta_grid = np.arcsin(rho_grid * np.sin(alpha))

        z_1_out = z_1(theta_grid[mask], phi_grid[mask])
        z_2_out = z_2(theta_grid[mask], phi_grid[mask])

        # difference in "waves"
        dz = (z_1_out - z_2_out)

        # wrap to [-0.5, 0.5) waves (shortest phase difference)
        dz_wrapped = ((dz + 0.5) % 1.0) - 0.5

        return np.sqrt(np.mean(dz_wrapped**2))

    #takes in the input from the booth correction algorithm
    def produce_correction_graphic(self, quality_data, new_corrections):
        a = self.microscope.alpha
        fig, axs = plt.subplots(2, 2, sharex = "col")
        
        many_composite(fig, axs[:, 0], self.raw_image_x, self.raw_image_y,
                       [self.get_uncorrected_img(), self.get_corrected_img()],
                       [self.get_empty_map(), self.get_corrected_map()], a,
                       is_edge = True, is_horizontal=False)

        #iterations
        iteration_no = np.arange(0, len(quality_data))
        #get quality data for the first image taken in each iteration
        qualities = quality_data[:, 0, 0]
        #get z_map for the correction in each iteration
        new_corrections = np.concatenate([[np.zeros_like(new_corrections[0])], new_corrections])
        z_maps = [create_zernike_function(self.innate_aberration.modes,
                                          self.innate_aberration.strengths + corrections, a) for corrections in new_corrections]
        errors = [SLM.zernike_RMS(z_map, self.get_empty_map(), a) for z_map in z_maps]
        
        axs[0][1].plot(iteration_no, qualities/self.ideal_quality, marker = "o")
        axs[0][1].set_ylabel("Image Quality")

        axs[1][1].plot(iteration_no, errors, marker = "o")
        axs[1][1].set_xlabel("Iteration No.")
        axs[1][1].set_ylabel("Correction RMS")
        axs[1][1].set_xticks(iteration_no)
