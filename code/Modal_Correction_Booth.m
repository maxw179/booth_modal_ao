clear;
clc;
close all;

%% === Booth 2002 Modal Correction ===

% Parameters
% n = 1.3225;                  % Refractive index (water)
% NA = 1.05;                  % Numerical aperture
% grid_size = 128;            % BFP resolution
% L_bfp = 15.12;              % BFP diameter [mm]
% 
% dx_bfp = L_bfp / grid_size;
% 
% [x_bfp, y_bfp] = meshgrid((-grid_size/2:grid_size/2-1) * dx_bfp);
% r_bfp = sqrt(x_bfp.^2 + y_bfp.^2);
% phi_bfp = atan2(y_bfp, x_bfp);
% 
% R_bfp = L_bfp / 2;
% theta = asin(r_bfp * NA / (n * R_bfp)); % Abbe condition
% theta(r_bfp > R_bfp) = 0;
% 
% x_Zernike = x_bfp / R_bfp;
% y_Zernike = y_bfp / R_bfp;
% rho = sqrt(x_Zernike.^2 + y_Zernike.^2); % Normalized r
% mask = rho <= 1; % Mask within aperture

% Arbitrary mode: Spherical
n = 2;
m = 0;
c = 1;

% Z = zeros(size(rho));
% Z(mask) = zernike(n, m, rho(mask), phi_bfp(mask));
% 
% W = c * Z;  % initial aberrated phase
% W = zeros(size(Z)); % non-aberrated wavefront

% --- Plot initial wavefront ---
% figure;
% imagesc(W); axis image; colorbar;
% title('Initial Wavefront W (Iteration 0)');
% drawnow;


I0 = PSF_Watzky(n, m, c); % FWHM of original image
                    % Using this measure temporarily, ideally it is signal

b = 0.5;                  % bias amplitude [radians]
gain = 0.6;               % correction gain (0.3–1)
n_iters = 10;              % number of correction iterations

a_est = -1; % initialize estimated amplitude

I_mat = []; % initialize intensity matrix
coeff_mat = []; % initialize coefficient matrix
I_mat(end+1) = I0;
coeff_mat(end+1) = a_est;

for iter = 1:n_iters
    % Compute positive and negative bias PSFs
    % W_pos = W + b * Z;
    % W_neg = W - b * Z;
    
    coeff_mat(end+1) = a_est + iter*0.001;
    coeff_mat(end+1) = a_est - iter*0.001;

    % Compute integrated FWHM (sub for signal insten.)
    I1 = PSF_Watzky(n, m, c + a_est  + iter*0.05);
    I2 = PSF_Watzky(n, m, c + a_est  - iter*0.05);
    
    I_mat(end+1) = I1;
    I_mat(end+1) = I2;

    dI = I1 - I2;               % Booth's delta_W signal
    % a_est = a_est - gain * dI;  % update estimate
    
    % % Apply correction by subtracting estimated aberration
    % W = W - gain * dI * Z;

    % --- Plot wavefront after each iteration ---
    % figure;
    % subplot(1,3,1);
    % imagesc(W .* mask); axis image; colorbar;
    % title(['W (Iteration ', num2str(iter), ')']);
    % 
    % subplot(1,3,2);
    % imagesc(W_pos .* mask); axis image; colorbar;
    % title('W + bZ (W_{pos})');
    % 
    % subplot(1,3,3);
    % imagesc(W_neg .* mask); axis image; colorbar;
    % title('W - bZ (W_{neg})');
    % 
    % drawnow;


    fprintf('delta_I = %.4f, new a_est = %.4f rad\n', dI, a_est);
end
fprintf('Final corrected amplitude: %.3f rad\n', a_est);

disp('coeff_mat = ');
disp(coeff_mat);

disp('I_mat = ');
disp(I_mat);

%% === Quadratic Fit to I(a) ==
% Fit I(a) = p1*a^2 + p2*a + p3
p = polyfit(coeff_mat, I_mat, 2);

% Theoretical maximum occurs at vertex of parabola:
a_max = -p(2)/(2*p(1));
I_max = polyval(p, a_max);

fprintf('\nQuadratic Fit Results:\n');
fprintf('a_max = %.4f rad\n', a_max);
fprintf('I_max = %.4f (normalized intensity)\n', I_max);

% Plot data and quadratic fit
a_fit = linspace(min(coeff_mat)-0.2, max(coeff_mat)+0.2, 200);
I_fit = polyval(p, a_fit);

figure;
plot(coeff_mat, I_mat, 'bo', 'MarkerFaceColor','b'); hold on;
plot(a_fit, I_fit, 'r-', 'LineWidth', 2);
plot(a_max, I_max, 'ks', 'MarkerFaceColor','y', 'MarkerSize', 10);
xlabel('Aberration coefficient (a) [rad]');
ylabel('Integrated intensity I(a)');
legend('Simulated data', 'Quadratic fit', 'Maximum');
title('Booth Modal Correction — Quadratic Fit to ΔI vs. a');
grid on;

I_mat;
coeff_mat;

%% Zernike function

function Z = zernike(n, m, rho, theta)
    % Check for valid input
    
    if abs(m) > n
        error('Invalid Zernike indices: n must be >= |m|');
    end

    % Compute the radial part
    R = zeros(size(rho));
    for k = 0:(n - abs(m))
        alpha = (-1)^k * factorial(2*n - abs(m) - k) / ...
            (factorial(k) * factorial(n - k) * factorial(n - abs(m) - k));
        R = R + alpha * rho.^(2*(n - k) - abs(m));
    end

    % Apply the angular part
    if m > 0
        Z = R .* cos(m * theta);
    elseif m < 0
        Z = R .* sin(-m * theta);  % sin(|m|*theta)
    else
        Z = R;
    end
end