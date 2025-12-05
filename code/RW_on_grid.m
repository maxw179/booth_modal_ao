function [PSF] = RW_on_grid(x, y, z, alpha, k, f, fwhm_pupil, n, nTheta, nPhi, aberrationFunc, N)
%
%   [PSF] = RW_on_grid(x, y, z, alpha, k, f, fwhm_pupil, n, ...
%                                 nTheta, nPhi, aberrationFunc, N)
%
%   All lengths (x,y,z,f,fwhm_pupil, wavelength used to define k) should
%   be in the SAME UNITS (ideally mm)
%
%   Inputs:
%       x, y, z       - VECTORIZED coordinates of field point (scalars, same units as f)
%       alpha         - maximum cone angle (radians), alpha = asin(NA / n)
%       k             - wavenumber in medium [1/unit], ideally k = 2*pi*n/lambda
%       f             - focal length of objective (same units as x,y,z)
%       fwhm_pupil    - Gaussian FWHM in back focal plane [same units as f]
%       n             - refractive index
%       nTheta        - number of theta samples
%       nPhi          - number of phi samples
%       aberrationFunc - function handle @(thetaGrid,phiGrid) -> phase (rad),
%                        or [] / omitted for no aberration
%       N              - microscope order (e.g. 3 for 3-photon)
%
%   Outputs:
%       Ex, Ey, Ez    - complex field components at (x,y,z)

    nx = numel(x);
    ny = numel(y);
    nz = numel(z);

    Ex = zeros(nx, ny, nz);
    Ey = zeros(nx, ny, nz);
    Ez = zeros(nx, ny, nz);

    fprintf('Computing PSF on %d x %d grid (y-plane only)...\n', nx, nz);

    for ix = 1:nx
        for iy = 1:ny   
            for iz = 1:nz
                [Ex(ix,iy,iz), Ey(ix,iy,iz), Ez(ix,iy,iz)] = ...
                    single_integrate_RW( ...
                        x(ix), y(iy), z(iz), ...
                        alpha, k, f, fwhm_pupil, n, ...
                        nTheta, nPhi, aberrationFunc);   % aberrationFunc = [] means no added aberration
            end
        end
    end

    PSF = abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2;
    PSF = PSF.^N;    % N-photon PSF
end

function [Ex, Ey, Ez] = single_integrate_RW(x, y, z, alpha, k, f, fwhm_pupil, n, nTheta, nPhi, aberrationFunc)
%
%   [Ex, Ey, Ez] = E_integrate_RW(x, y, z, alpha, k, f, fwhm_pupil, n, ...
%                                 nTheta, nPhi, aberrationFunc)
%
%   All lengths (x,y,z,f,fwhm_pupil, wavelength used to define k) should
%   be in the SAME UNITS (ideally mm)
%
%   Inputs:
%       x, y, z       - single coordinates of field point (scalars, same units as f)
%       alpha         - maximum cone angle (radians), alpha = asin(NA / n)
%       k             - wavenumber in medium [1/unit], ideally k = 2*pi*n/lambda
%       f             - focal length of objective (same units as x,y,z)
%       fwhm_pupil    - Gaussian FWHM in back focal plane [same units as f]
%       n             - refractive index
%       nTheta        - number of theta samples
%       nPhi          - number of phi samples
%       aberrationFunc - function handle @(thetaGrid,phiGrid) -> phase (rad),
%                        or [] / omitted for no aberration
%
%   Outputs:
%       Ex, Ey, Ez    - complex field components at (x,y,z)

    %defualts if nPhi, nTheta, or aberationFunc are not provided
    if nargin < 11 || isempty(nPhi)
        nPhi = 128;
    end
    if nargin < 10 || isempty(nTheta)
        nTheta = 128;
    end
    if nargin < 12
        aberrationFunc = [];
    end

    theta = linspace(0, alpha, nTheta);
    phi   = linspace(0, 2*pi, nPhi);
    dTheta = theta(2) - theta(1);
    dPhi   = phi(2)   - phi(1);
    [phiGrid, thetaGrid] = meshgrid(phi, theta);

    cosT = cos(thetaGrid);
    sinT = sin(thetaGrid);
    amp = gaussian_amplitude(thetaGrid, f, n, fwhm_pupil); 
    inside_factor = amp .* sqrt(cosT) .* sinT;

    %handling added aberration
    if isempty(aberrationFunc)
        aberrationPhase = 0;
    else
        aberrationPhase = aberrationFunc(thetaGrid, phiGrid);
    end

    defaultPhase = aberration_free_phase(k, x, y, z, thetaGrid, phiGrid);
    phase = exp(1i * (defaultPhase + aberrationPhase));

    [ax, ay, az] = strength_angular(thetaGrid, phiGrid);

    pre_factor = -1i * k * f / (2*pi);

    integrand_x = inside_factor .* ax .* phase;
    integrand_y = inside_factor .* ay .* phase;
    integrand_z = inside_factor .* az .* phase;

    Ex = pre_factor * sum(integrand_x, 'all') * dTheta * dPhi;
    Ey = pre_factor * sum(integrand_y, 'all') * dTheta * dPhi;
    Ez = pre_factor * sum(integrand_z, 'all') * dTheta * dPhi;
end

function amp = gaussian_amplitude(theta, f, n, fwhm_pupil)
    r_bfp = f .* n .* sin(theta);  % mm if f is mm
    amp = exp(-4*log(2) .* (r_bfp.^2) ./ (fwhm_pupil.^2));
end

function phi_tot = aberration_free_phase(k, x, y, z, theta, phi)
    phi_tot = k .* ( ...
        x .* sin(theta) .* cos(phi) + ...
        y .* sin(theta) .* sin(phi) + ...
        z .* cos(theta) );
end


function [a_x, a_y, a_z] = strength_angular(theta, phi)
    cosT = cos(theta);
    sinT = sin(theta);
    cosP = cos(phi);
    sinP = sin(phi);

    a_x = cosT + (1 - cosT) .* (sinP.^2);
    a_y = (cosT - 1) .* cosP .* sinP;
    a_z = -sinT .* cosP;
end
