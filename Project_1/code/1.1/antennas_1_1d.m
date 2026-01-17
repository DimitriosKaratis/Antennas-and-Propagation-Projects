%% Antennas and Propagation Project - Nov 2025
% KARATIS DIMITRIOS 10775 
%% 1.1d


% Parameters
clc; close all; clear;   % Clear workspace and figures
f = 1e9;                 % Frequency in Hz
c0 = 3e8;                % Speed of light (m/s)
lambda = c0 / f;         % Wavelength in meters
d = lambda / 2;          % Distance between dipoles
Nx = 24;                 % Number of elements along x-axis
Nz = 12;                 % Number of elements along z-axis
k = 2 * pi / lambda;     % Wave number

% Dipole constants
I0 = 1;                  % Current amplitude in Ampere
h0 = 120*pi;             % Impedance of free space
r = 50;                  % Distance to observation point in meters
l = lambda / 2;          % Dipole length

% Max radiation angles (theta, phi)
theta_m = [pi/2, pi/3];
phi_m = [pi/2, pi/3, pi/6];


% Define mesh grid for theta and phi (180 points for theta, 360 for phi)
Theta_linspace = linspace(0, pi, 180);
Phi_linspace = linspace(0, 2*pi, 360);

[theta, phi] = meshgrid(Theta_linspace, Phi_linspace);

dtheta = Theta_linspace(2)-Theta_linspace(1);
dphi   = Phi_linspace(2)-Phi_linspace(1);


% Initialize the electric field matrix
E_abs = zeros(size(phi));

% Output table
results = [];

fprintf('1.1d\n')
% Loop for each combination of theta_m and phi_m
for i = 1:length(theta_m)
    for j = 1:length(phi_m)
        
        % Calculate phase shifts for x-axis and z-axis
        delta_x = calculate_delta_x(phi_m(j), theta_m(i), k, d, 0);
        delta_z = calculate_delta_z(theta_m(i), k, d, 0);

        % Electric field calculation for horizontal (phi) and vertical (theta)
        E0_abs = abs(h0 * I0 * k * l * sin(theta) / (4*pi*r)); % Electric field constant
        
        % Calculate the phase shifts in x and z directions
        psi_x = k * d * cos(phi) .* sin(theta) + delta_x;
        psi_z = k * d * cos(theta) + delta_z;
        
        % Calculate the Array Factor (AF)
        AF_abs = abs(sin(Nx * psi_x / 2) ./ sin(psi_x / 2)) .* abs(sin(Nz * psi_z / 2) ./ sin(psi_z / 2));
        
        % Total electric field
        E_abs = AF_abs .* E0_abs;   

        % Radiation intensity
        U = E_abs.^2;
    
        % Numerical integration for total radiated power (Riemann Sum)
        Prad = sum(sum( U .* sin(theta) )) * dtheta * dphi;
    
        % Max radiation
        Umax = max(U(:));
    
        % Directivity
        D = 4*pi * Umax / Prad;
        D_dBi = 10*log10(D);

        % Save results
        results = [results;
                   theta_m(i)*180/pi, phi_m(j)*180/pi, D, D_dBi];    
        
    end
end

% Display results
fprintf('---------------------------------------------\n');
disp('   θ_{max}   φ_{max}  Direct.    Direct.(dBi)');
disp(results);



%% Functions to Calculate Phase Shifts
% Calculate delta for x-axis (phase shift along x-direction)
function delta_x = calculate_delta_x(phi, theta, k, dx, psi_x)
    delta_x = psi_x - k * dx * cos(phi) * sin(theta); 
end

% Calculate delta for z-axis (phase shift along z-direction)
function delta_z = calculate_delta_z(theta, k, dz, psi_z)
    delta_z = psi_z - k * dz * cos(theta); 
end

