%% Antennas and Propagation Project - Nov 2025
% KARATIS DIMITRIOS 10775 
%% 1.1a

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

% Horizontal radiation pattern figure
figure;

% Set up subplots
subplot_row = length(theta_m);
subplot_col = length(phi_m);

% Loop for horizontal pattern
for i = 1:length(theta_m)
    for j = 1:length(phi_m)
        % Compute phase shifts
        delta_x = calculate_delta_x(phi_m(j), theta_m(i), k, d, 0);
        delta_z = calculate_delta_z(theta_m(i), k, d, 0);

        % Small delta correction
        if abs(delta_x) < 1e-7
            delta_x = 0;
        end
        if abs(delta_z) < 1e-7
            delta_z = 0;
        end

        % Horizontal pattern calculations
        phi_hor = linspace(0, 2*pi, 360);
        theta_hor = pi/2;

        E0_abs_hor = abs(h0 * I0 * k * l * sin(theta_hor) / (4*pi*r));

        psi_x = k * d * cos(phi_hor) * sin(theta_hor) + delta_x;
        psi_z = k * d * cos(theta_hor) + delta_z;

        AF_abs_hor = abs(sin(Nx * psi_x / 2) ./ sin(psi_x / 2)) .* abs(sin(Nz * psi_z / 2) ./ sin(psi_z / 2));
        
        E_abs_hor = AF_abs_hor .* E0_abs_hor;
        E_abs_hor_norm = E_abs_hor / max(E_abs_hor);


        % Plot horizontal radiation
        subplot(subplot_row, subplot_col, (i-1)*subplot_col + j);
        polarplot(phi_hor, E_abs_hor_norm);
        title(['Horizontal Pattern (θ_{max} = ', num2str(theta_m(i)*(180/pi)), '°, φ_{max} = ', num2str(phi_m(j)*(180)/pi), '°)']);
    end
end



% Vertical radiation pattern figure
figure;

% Loop for vertical pattern
for i = 1:length(theta_m)
    for j = 1:length(phi_m)
        % Compute phase shifts
        delta_x = calculate_delta_x(phi_m(j), theta_m(i), k, d, 0);
        delta_z = calculate_delta_z(theta_m(i), k, d, 0);

        % Small delta correction
        if abs(delta_x) < 1e-7
            delta_x = 0;
        end
        if abs(delta_z) < 1e-7
            delta_z = 0;
        end

        % Vertical pattern calculations
        theta_ver = linspace(0, 2*pi, 360);
        phi_ver = phi_m(j);
        E0_abs_ver = abs(h0 * I0 * k * l * sin(theta_ver) / (4*pi*r));

        psi_x_ver = k * d * cos(phi_ver) * sin(theta_ver) + delta_x;
        psi_z_ver = k * d * cos(theta_ver) + delta_z;

        AF_abs_ver = abs(sin(Nx * psi_x_ver / 2) ./ sin(psi_x_ver / 2)) .* abs(sin(Nz * psi_z_ver / 2) ./ sin(psi_z_ver / 2));
        
        E_abs_ver = AF_abs_ver .* E0_abs_ver;
        E_abs_ver_norm = E_abs_ver / max(E_abs_ver);

        % Plot vertical radiation
        subplot(subplot_row, subplot_col, (i-1)*subplot_col + j);
        polarplot(theta_ver, E_abs_ver_norm);
        title(['Vertical Pattern (θ_{max} = ', num2str(theta_m(i)*(180/pi)), '°, φ_{max} = ', num2str(phi_m(j)*(180)/pi), '°)']);
    end
end


%% Functions to Calculate Phase Shifts
% Calculate delta for x-axis
function delta_x = calculate_delta_x(phi, theta, k, dx, psi_x)
    delta_x = psi_x - k * dx * cos(phi) * sin(theta);           
end

% Calculate delta for z-axis
function delta_z = calculate_delta_z(theta, k, dz, psi_z)
    delta_z = psi_z - k * dz * cos(theta); 
end


%% EXTRA WORK: PLOTS FOR VALIDATION PURPOSES (VISUALLY CHECK THAT |E| = |E0|*|AFxz|, which is correct!)

% Figure comparing E0 for horizontal and vertical patterns
figure;
E0_abs_ver = abs(h0 * I0 * k * l * sin(theta_ver) / (4*pi*r));
E0_abs_ver_norm = E0_abs_ver / max(E0_abs_ver);
E0_abs_hor = abs(h0 * I0 * k * l / (4*pi*r));
E0_abs_hor_norm = E0_abs_hor / max(E0_abs_hor);

% Plot E0 for vertical and horizontal dipoles
subplot(1,2,1);
polarplot(theta_ver, E0_abs_ver_norm);
title('Vertical Pattern for E_0');
subplot(1,2,2);
phi_hor = linspace(0, 2*pi, 360);
polarplot(phi_hor, E0_abs_hor_norm * ones(size(phi_hor)));
title('Horizontal Pattern for E_0');

% Adjust layout for better appearance
sgtitle('Dipole Radiation Patterns');




% Final validation of Array Factor (AF) and Electric Field
figure;
for i = 1:length(theta_m)
    for j = 1:length(phi_m)
        % Compute phase shifts
        delta_x = calculate_delta_x(phi_m(j), theta_m(i), k, d, 0);
        delta_z = calculate_delta_z(theta_m(i), k, d, 0);

        % Small delta correction
        if abs(delta_x) < 1e-7
            delta_x = 0;
        end
        if abs(delta_z) < 1e-7
            delta_z = 0;
        end

        % Vertical pattern validation
        theta_ver = linspace(0, 2*pi, 360);
        phi_ver = phi_m(j);

        E0_abs_ver = abs(h0 * I0 * k * l * sin(theta_ver) / (4*pi*r));

        psi_x_ver = k * d * cos(phi_ver) * sin(theta_ver) + delta_x;
        psi_z_ver = k * d * cos(theta_ver) + delta_z;

        AF_abs_ver = abs(sin(Nx * psi_x_ver / 2) ./ sin(psi_x_ver / 2)) .* abs(sin(Nz * psi_z_ver / 2) ./ sin(psi_z_ver / 2));
        AF_abs_ver_norm = AF_abs_ver / max(AF_abs_ver);

        % Plot AF validation
        subplot(subplot_row, subplot_col, (i-1)*subplot_col + j);
        polarplot(theta_ver, AF_abs_ver_norm);
        title(['Vertical Pattern (θ_{max} = ', num2str(theta_m(i)*(180/pi)), '°, φ_{max} = ', num2str(phi_m(j)*(180)/pi), '°)']);
    end
end
sgtitle('Polar plots of |AF_{xz}| for validation');
