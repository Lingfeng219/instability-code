function [B, A, om, phi] = baroclinic_eigen_mode(u, v, rho, k, l, lat, z, alphax, alphay, rho0, visc, qx, qy)
% baroclinic_eigen_mode
% Solves the linear normal mode instability problem considering submesoscale 
% eddy mixing (in terms of viscosity, ref. Vollmer 2013), using the matrix 
% method (Smith et al., 2007).
%
% Usage:
% [B, A, om, phi] = baroclinic_eigen_mode(u, v, rho, k, l, lat, z, alphax, alphay, rho0, visc, qx, qy)
%
% Inputs:
%   u, v   - Velocity component profiles
%   rho    - Potential density profile
%   k, l   - Wavenumbers in x and y directions
%   lat    - Latitude (degrees)
%   z      - Vertical coordinate vector (evenly spaced)
%   alphax, alphay - Bottom topography slopes in x and y directions
%   rho0   - Reference density (default: 1030 kg/m³)
%   visc   - Horizontal viscosity profile or constant
%   qx, qy - Horizontal vorticity gradient (default: 0)
%
% Outputs:
%   B, A   - Matrices for the generalized eigenvalue problem: om * B * Phi = A * Phi
%   om     - Complex growth rate
%   phi    - Complex stream function eigenfunction
%
% Reference:
% Feng 2024, JGR, Global Distribution and Seasonal Variations of 
% Charney-type Submesoscale Baroclinic Instabilities (C-SBCIs)
%
% Author: Chuanyu Liu, IOCAS
% Date:   2025/03/19 (Updated for GitHub release)

% Set default values for missing inputs
if nargin < 13, qx = 0; qy = 0; end
if nargin < 11, visc = 0; end
if nargin < 10, rho0 = 1030; end
if nargin < 8, alphax = 0; alphay = 0; end

% Constants
small = 1e-20;
R = 6.378e6; % Earth's radius (m)
f = 1.458e-4 * sind(lat); % Coriolis parameter
beta = 1.458e-4 / R * cosd(lat); % Rossby parameter (df/dy)
grav = 9.81; % Gravitational acceleration (m/s²)
K2 = k^2 + l^2; % Square of total wavenumber

% Compute vertical grid spacing (dz)
n = length(z);
dz = diff(z);
dz(n) = dz(n - 1); % Assign last dz same as previous

% Ensure unique density profile
[B_unique, ~] = unique(rho);
for bb = 1:length(B_unique)
    idb = find(rho == B_unique(bb));
    if length(idb) > 1
        for bb_idx = 2:length(idb)
            rho(idb(bb_idx)) = rho(idb(1)) + 1e-3 * bb_idx;
        end
    end
end

% Construct T-matrix (from Smith 2007, Eq. 2.1 & B.1)
T = zeros(n, n);
for m = 2:n-1
    T(m, m-1) = 1 / ((rho(m) - rho(m-1)) * dz(m));
    T(m, m)   = - (T(m, m-1) + 1 / ((rho(m+1) - rho(m)) * dz(m)));
    T(m, m+1) = 1 / ((rho(m+1) - rho(m)) * dz(m));
end

% Boundary conditions for T
T(1, 1) = -1 / ((rho(2) - rho(1)) * dz(1));
T(1, 2) = 1 / ((rho(2) - rho(1)) * dz(1));
T(n, n-1) = 1 / ((rho(n) - rho(n-1)) * dz(n));
T(n, n) = -1 / ((rho(n) - rho(n-1)) * dz(n));
T = T * (f^2 * rho0 / grav);

% Construct B matrix
B = T;
for m = 1:n
    B(m, m) = T(m, m) - K2;
end

% Compute Qx and Qy
qx(isnan(qx)) = 0;
qy(isnan(qy)) = 0;
Qx = T * v + qx;
Qy = beta - T * u - qy;

% Construct diagonal matrices
Qxm = diag(Qx);
Qym = diag(Qy);

% Incorporate bottom topography slope effects
Qxm(n, n) = Qxm(n, n) + f * alphax / dz(n);
Qym(n, n) = Qym(n, n) + f * alphay / dz(n);

% Compute A2 matrix
A2 = k * Qym - l * Qxm;

% Construct velocity and viscosity matrices
Um = diag(u);
Vm = diag(v);

if length(visc) == n
    Ahm = diag(visc);
else
    Ahm = diag(visc * ones(n, 1));
end

% Compute AL matrix
AL = k * Um + l * Vm - 1i * Ahm * K2;

% Compute A1 and final A matrix
A1 = AL * B;
A = A2 + A1;

% Solve the generalized eigenvalue problem
try 
    eig(A, B);
    [phi, om] = eig(A, B);
catch
    phi = [];
    om = [];
end

% End of function
end
