%% run_RW_PSF_port.m
% Wrapper to compute PSF using E_integrate_RW

clear; clc; close all;

%% Inputs / parameters

N = 3;                      % N-photon order
fwhm_pupil_mm = 8.24;       % mm, Gaussian FWHM at BFP

wavelength_um = 1.3;        % µm 
wavelength_mm = wavelength_um / 1000;  % mm

% x,y,z coordinates relative to focus, in µm
x_um = -2:0.05:2;           % µm
y_um = 0;                   % single line
z_um = -3:0.05:3;           % µm

% Objective parameters
tubeLength_mm  = 200;       % mm
NA             = 1.05;
magnification  = 25;
n_medium       = 1.3225;    % refractive index in sample

% Derived optical parameters
f_mm   = tubeLength_mm / magnification;        % focal length in mm
alpha  = asin(NA / n_medium);                  % max cone angle
k      = 2*pi * n_medium / wavelength_mm;      % wavenumber in medium [1/mm]

% Debye integration sampling
nTheta = 128;
nPhi   = 128;

%% Convert coordinates to mm for integra

x_mm = x_um / 1000;     % µm -> mm
y_mm = y_um / 1000;
z_mm = z_um / 1000;


%% Compute PSF
PSF = RW_on_grid(x_mm, y_mm, z_mm, alpha, k, f_mm, fwhm_pupil_mm, n_medium, nTheta, nPhi, [], N); %no aberration function

%% Plot PSF: axial and lateral

% Axial PSF at x closest to 0
[~, ix0] = min(abs(x_um - 0));
profile_z = PSF(ix0, :);
profile_z = profile_z / max(profile_z);

figure;
plot(z_um, profile_z, 'LineWidth', 1.5);
xlabel('z (\mum)');
ylabel(sprintf('%d-photon PSF (normalized)', N));
title('Axial PSF (x \approx 0, y = 0)');
grid on;

% Lateral PSF at z closest to 0
[~, iz0] = min(abs(z_um - 0));
profile_xy = PSF(:, iz0);
profile_xy = profile_xy / max(profile_xy);

figure;
plot(x_um, profile_xy, 'LineWidth', 1.5);
xlabel('x (\mum)');
ylabel(sprintf('%d-photon PSF (normalized)', N));
title('Lateral PSF (z \approx 0, y = 0)');
grid on;

%% FWHM calculations 

% Axial FWHM (z)
flagz = profile_z >= 0.5;
FWHMz = sum(flagz) * (z_um(2) - z_um(1));   % µm

% Lateral FWHM (x)
flagxy = profile_xy >= 0.5;
FWHMxy = sum(flagxy) * (x_um(2) - x_um(1)); % µm

% Optional bead deconvolution
FWHMbead = 0.2;  % µm
FWHMxy_meas = sqrt(FWHMxy.^2 + FWHMbead.^2);
FWHMz_meas  = sqrt(FWHMz.^2  + FWHMbead.^2);

fprintf('FWHM_xy  = %.3f µm (raw),  %.3f µm (with 0.2 µm bead)\n', FWHMxy, FWHMxy_meas);
fprintf('FWHM_z   = %.3f µm (raw),  %.3f µm (with 0.2 µm bead)\n', FWHMz, FWHMz_meas);
