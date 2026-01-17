%% Antennas and Propagation Project - Nov 2025
% KARATIS DIMITRIOS 10775 
%% 1.1b

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
[theta, phi] = meshgrid(linspace(0, pi, 180), linspace(0, 2*pi, 360));



% Initialize the electric field matrix
E_abs = zeros(size(phi));

% Loop for each combination of theta_m and phi_m
for i = 1:length(theta_m)
    for j = 1:length(phi_m)
        
        % Calculate phase shifts for x-axis and z-axis
        delta_x = calculate_delta_x(phi_m(j), theta_m(i), k, d, 0);
        delta_z = calculate_delta_z(theta_m(i), k, d, 0);

        % Small delta correction
        if abs(delta_x) < 1e-7
            delta_x = 0;
        end
        if abs(delta_z) < 1e-7
            delta_z = 0;
        end
        
        % Electric field calculation for horizontal (phi) and vertical (theta)
        E0_abs = abs(h0 * I0 * k * l * sin(theta) / (4*pi*r)); % Electric field constant
        
        % Calculate the phase shifts in x and z directions
        psi_x = k * d * cos(phi) .* sin(theta) + delta_x;
        psi_z = k * d * cos(theta) + delta_z;
        
        % Calculate the Array Factor (AF)
        AF_abs = abs(sin(Nx * psi_x / 2) ./ sin(psi_x / 2)) .* abs(sin(Nz * psi_z / 2) ./ sin(psi_z / 2));
        
        % Total electric field
        E_abs = AF_abs .* E0_abs;
        E_abs_norm = E_abs / max(E_abs(:)); % Normalize the electric field

        % Convert spherical coordinates to Cartesian coordinates for plotting
        X = E_abs_norm .* sin(theta) .* cos(phi);  % X = r * sin(theta) * cos(phi)
        Y = E_abs_norm .* sin(theta) .* sin(phi);  % Y = r * sin(theta) * sin(phi)
        Z = E_abs_norm .* cos(theta);              % Z = r * cos(theta)
        
        % Plot the radiation pattern as a 3D surface
        figure;
        surf(X, Y, Z, E_abs_norm, 'EdgeColor', 'none');  % E_abs_norm for color
        colormap jet;                                    % Use jet colormap for better visualization
        colorbar;                       
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        title('');
        title(['3D Polar Radiation Pattern for (θ_{max} = ', num2str(theta_m(i)*(180/pi)), '°, φ_{max} = ', num2str(phi_m(j)*(180)/pi), '°)']);
        axis equal;    
        view(3);        

    end
end




%% Functions to Calculate Phase Shifts
% Calculate delta for x-axis (phase shift along x-direction)
function delta_x = calculate_delta_x(phi, theta, k, dx, psi_x)
    delta_x = psi_x - k * dx * cos(phi) * sin(theta); 
end

% Calculate delta for z-axis (phase shift along z-direction)
function delta_z = calculate_delta_z(theta, k, dz, psi_z)
    delta_z = psi_z - k * dz * cos(theta); 
end