%% Antennas and Propagation Project - Nov 2025
% KARATIS DIMITRIOS 10775 
%% 1.1st

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
theta_m = pi/2;
phi_m = 0;

disp('1.1st');
% New radiation angles (after the coordinate transformation)
theta_0_tonos = new_theta(theta_m, phi_m);
phi_0_tonos = new_phi(theta_m, phi_m);
print_old_and_new_angles(theta_m, phi_m, theta_0_tonos, phi_0_tonos);

% Define mesh grid for theta and phi (180 points for theta, 360 for phi)
Theta_linspace = linspace(0, pi, 180);
Phi_linspace = linspace(0, 2*pi, 360);

[theta, phi] = meshgrid(Theta_linspace, Phi_linspace);

dtheta = Theta_linspace(2)-Theta_linspace(1);
dphi   = Phi_linspace(2)-Phi_linspace(1);


% Initialize the electric field matrix
E_abs = zeros(size(phi));

% Loop for each combination of theta_m and phi_m
for i = 1:length(theta_m)
    for j = 1:length(phi_m)
        
        % Calculate phase shifts for x-axis and z-axis
        delta_z = calculate_delta_z(theta_m(i), k, d, 0);
        delta_x = calculate_delta_Hansen_Woodyard(k, d, Nx); % Hansen-Woodyard in x(y')
        
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
        
        %% 3D POLAR PLOT
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

        %% THEORETICAL DIRECTIVITY CALCULATION
        % Radiation intensity
        U = E_abs.^2;
    
        % NNumerical integration for total radiated power (Riemann Sum)
        Prad = sum(sum( U .* sin(theta) )) * dtheta * dphi;
    
        % Max radiation
        Umax = max(U(:));
    
        % Directivity
        D = 4*pi * Umax / Prad;
        D_dBi = 10*log10(D);


        fprintf("Method 1 (Theoretical) | The directivity D for (θ_{max}, φ_{max}) = (%.1f, %.1f) is: %.3f dBi\n", ...
                 theta_m(i)*180/pi, phi_m(j)*180/pi, D_dBi);

        %% DIRECTIVITY CALCULATION USING D
        % Cannot be used, since cos(θ') = 0
        
        fprintf('Method 2 (using Dx,Dz) | CANNOT BE USED\n');
        
        %% DIRECTIVITY CALCULATION USING HPBW
        
        Theta_z = calculate_HPBW_for_broadside(lambda, Nz, d); % (eyripleri ston z (x'))
        Theta_x = calculate_HPBW_for_Hansen_Woodyard(lambda, Nx, d); % Hansen-Woodyard ston x (y') (palios x, neos y diladi)
          
        D_hpbw = calculate_directivity_from_HPBW(Theta_z, Theta_x);
        D_hpbw_dbi = 10*log10(D_hpbw);

        fprintf('Method 3 (using HPBW)  | The directivity D for (θ_{max}, φ_{max}) = (%.1f, %.1f) is: %.3f dBi\n', ...
        theta_m(i)*(180/pi), phi_m(j)*(180/pi), D_hpbw_dbi);
    end
end







%% Helper Functions 

% Calculate delta for z-axis (phase shift along z-direction)
function delta_z = calculate_delta_z(theta, k, dz, psi_z)
    delta_z = psi_z - k * dz * cos(theta); 
end

% Calculate delta for z-axis (phase shift along z-direction)
function delta = calculate_delta_Hansen_Woodyard(k, d, N)
    delta = - k * d - (2.92/N); 
end

% This function calculates the new elevation angle (theta')
function theta_prime = new_theta(theta, phi)
    theta_prime = acos(sin(theta) * sin(phi));  % θ0' = cos⁻¹(sin(θ) * sin(φ))
end

% This function calculates the new azimuth angle (phi')
function phi_prime = new_phi(theta, phi)
    % In case we have 0/0
    if sin(theta) * cos(phi) < 10^(-4) && cos(theta) < 10^(-4) 
        phi_prime = pi/2; 
    else
        phi_prime = atan2(sin(theta) * cos(phi), cos(theta));  % φ0' = atan2(sin(θ) * cos(φ), cos(θ))
    end
end

% This function calculates Half Power Beamwidth (HPBW) for a broadside array
function HPBW = calculate_HPBW_for_broadside(lambda, N, d)
     HPBW = 48.4 * lambda / (N * d);
end

% This function calculates Half Power Beamwidth (HPBW) for the Hansen-Woodyard case
function HPBW = calculate_HPBW_for_Hansen_Woodyard(lambda, N, d)
     HPBW = 2 * rad2deg(acos(1 - 0.1398 * (lambda / (N*d))));
     
end

% This function calculates the directivity from HPBW values for an end-fire
function D_HPBW = calculate_directivity_from_HPBW(Theta_x_tonos, Theta_y_tonos)  
    D_HPBW = 32400 / (Theta_x_tonos * Theta_y_tonos);
end

% Function to print the old and new theta and phi
function print_old_and_new_angles(theta_m, phi_m, theta_0_tonos, phi_0_tonos)
    
    fprintf('----------------------------------------------------------------------\n');
    
    % Loop through each combination of theta and phi
    for i = 1:length(theta_m)
        for j = 1:length(phi_m)
            % Old angles
            theta_old = theta_m(i);
            phi_old = phi_m(j);
            
            % New angles from the calculated values
            theta_new = theta_0_tonos(i, j);
            phi_new = phi_0_tonos(i, j);
                      
            % Print the results
            fprintf('Old (θ, φ) = (%.2f, %.2f)° ---> ', ...
                    rad2deg(theta_old), rad2deg(phi_old));
            fprintf('New (θ, φ) = (%.2f, %.2f)°\n\n', ...
                    rad2deg(theta_new), rad2deg(phi_new));

        end
    end
end
