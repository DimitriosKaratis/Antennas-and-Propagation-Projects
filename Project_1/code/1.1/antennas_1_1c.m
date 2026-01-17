%% Antennas and Propagation Project - Nov 2025
% KARATIS DIMITRIOS 10775 
%% 1.1c

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

disp('1.1c');
% New radiation angles (after the coordinate transformation)
theta_0_tonos = [new_theta(pi/2, pi/2) new_theta(pi/2, pi/3) new_theta(pi/2, pi/6); 
                 new_theta(pi/3, pi/2) new_theta(pi/3, pi/3) new_theta(pi/3, pi/6)];

phi_0_tonos = [new_phi(pi/2, pi/2) new_phi(pi/2, pi/3) new_phi(pi/2, pi/6); 
               new_phi(pi/3, pi/2) new_phi(pi/3, pi/3) new_phi(pi/3, pi/6)];

print_old_and_new_angles(theta_m, phi_m, theta_0_tonos, phi_0_tonos);
fprintf('----------------------------------------------------------------------\n');
fprintf('Directivity Results\n\n');

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
            
       
        % For old z, new x' axis
        HPBW_broadside_z = 9;   % degrees (from Array Length diagram, Arr Len= Nz*dz/λ=6 and 90 deg from end-fire)
        HPBW60_z  = 10.5;       % degrees (from Array Length diagram, Arr Len= Nz*dz/λ=6 and 60 deg from end-fire)
        HPBW40_9_z  = 14.0;     % degrees (from Array Length diagram, Arr Len= Nz*dz/λ=6 and 40.89 deg from end-fire)
        HPBW56_3_z  = 11.5;     % degrees (from Array Length diagram, Arr Len= Nz*dz/λ=6 and 56.31 deg from end-fire)
        
        % For old x, new y' axis
        HPBW_broadside_x = 4.8;  % degrees (from Array Length diagram, Arr Len= Nx*dx/λ=12 and 90 deg from end-fire)
        HPBW60_x = 5.4;          % degrees (from Array Length diagram, Arr Len= Nx*dx/λ=12 and 60 deg from end-fire)
        HPBW30_x = 9.5;          % degrees (from Array Length diagram, Arr Len= Nx*dx/λ=12 and 30 deg from end-fire)
        HPBW48_6_x = 6.7;        % degrees (from Array Length diagram, Arr Len= Nx*dx/λ=12 and 90-41.41 = 48.6 deg from end-fire)
        HPBW25_7_x = 12.5;       % degrees (from Array Length diagram, Arr Len= Nx*dx/λ=12 and 90-64.34 = 25.66 deg from end-fire)


        if theta_m(i) == pi/2
           Dz_tonos = calculate_directivity_from_D_for_broadside(lambda, Nz, d);     % Broadside in z(x')

           if phi_m(j) == pi/2
               Dx_tonos = calculate_directivity_from_D_for_broadside(lambda, Nx, d); % Broadside in x(y')

           elseif phi_m(j) == pi/3 
               Dx = calculate_directivity_from_D_for_broadside(lambda, Nx, d); 
               Theta_broad_x = HPBW_broadside_x;
               Theta_x = HPBW60_x; 
               Dx_tonos = Dx * Theta_broad_x / Theta_x;

           else
               Dx = calculate_directivity_from_D_for_broadside(lambda, Nx, d); 
               Theta_broad_x = HPBW_broadside_x;
               Theta_x = HPBW30_x; 
               Dx_tonos = Dx * Theta_broad_x / Theta_x;

           end

        elseif theta_m(i) == pi/3
           if phi_m(j) == pi/2
               Dx_tonos = calculate_directivity_from_D_for_broadside(lambda, Nx, d); % Broadside in x(y')
               
               Dz = calculate_directivity_from_D_for_broadside(lambda, Nz, d); 
               Theta_broad_z = HPBW_broadside_z;
               Theta_z = HPBW60_z; 
               Dz_tonos = Dz * Theta_broad_z / Theta_z;

           elseif phi_m(j) == pi/3
               Dx = calculate_directivity_from_D_for_broadside(lambda, Nx, d); 
               Theta_broad_x = HPBW_broadside_x;
               Theta_x = HPBW48_6_x; 
               Dx_tonos = Dx * Theta_broad_x / Theta_x;

               Dz = calculate_directivity_from_D_for_broadside(lambda, Nz, d); 
               Theta_broad_z = HPBW_broadside_z;
               Theta_z = HPBW40_9_z; 
               Dz_tonos = Dz * Theta_broad_z / Theta_z; 

           else
               Dx = calculate_directivity_from_D_for_broadside(lambda, Nx, d); 
               Theta_broad_x = HPBW_broadside_x;
               Theta_x = HPBW25_7_x; 
               Dx_tonos = Dx * Theta_broad_x / Theta_x;

               Dz = calculate_directivity_from_D_for_broadside(lambda, Nz, d); 
               Theta_broad_z = HPBW_broadside_z;
               Theta_z = HPBW56_3_z; 
               Dz_tonos = Dz * Theta_broad_z / Theta_z;                 
           end
        end
        
        D_total = pi* cos(theta_0_tonos(i,j)) * Dx_tonos * Dz_tonos;
        D_total = 10*log10(D_total);

        fprintf('Method 1 (using Dx,Dz) | The directivity D for (θ_{max}, φ_{max}) = (%.1f, %.1f) is: %.3f dBi\n', ...
        theta_m(i)*(180/pi), phi_m(j)*(180/pi), D_total);
        
    end
end
fprintf('------------------------------------------------------------------------------------------------\n');





% Loop for each combination of theta_m and phi_m
for i = 1:length(theta_m)
    for j = 1:length(phi_m)
        
        % Old angles
        theta_old = theta_m(i);
        phi_old = phi_m(j);
                
        % New angles from the calculated values
        theta_new = theta_0_tonos(i, j);
        phi_new = phi_0_tonos(i, j);
                
        if theta_m(i) == pi/2
           Theta_z = HPBW_broadside_z;        % Broadside in z(x')
           
           if phi_m(j) == pi/2
                Theta_x = HPBW_broadside_x;   % Broadside in x(y')
               
           elseif phi_m(j) == pi/3 
                Theta_x = HPBW60_x;  
           else
                Theta_x = HPBW30_x;  
           end

        elseif theta_m(i) == pi/3

           if phi_m(j) == pi/2
               Theta_x = HPBW_broadside_x;   % Broadside in x(y')
               Theta_z = HPBW60_z;   

           elseif phi_m(j) == pi/3
               Theta_x = HPBW48_6_x; 
               Theta_z = HPBW40_9_z;
           else
               Theta_x = HPBW25_7_x; 
               Theta_z = HPBW56_3_z; 
           end
        end


        D_hpbw = calculate_directivity_from_HPBW(theta_new, phi_new, Theta_z, Theta_x);
        D_hpbw_dbi = 10*log10(D_hpbw);

        fprintf('Method 2 (using HPBW) | The directivity D for (θ_{max}, φ_{max}) = (%.1f, %.1f) is: %.3f dBi\n', ...
        theta_m(i)*(180/pi), phi_m(j)*(180/pi), D_hpbw_dbi);
    end
end




%% Helper Functions 

% Calculate delta for x-axis (phase shift along x-direction)
function delta_x = calculate_delta_x(phi, theta, k, dx, psi_x)
    delta_x = psi_x - k * dx * cos(phi) * sin(theta); 
end


% Calculate delta for z-axis (phase shift along z-direction)
function delta_z = calculate_delta_z(theta, k, dz, psi_z)
    delta_z = psi_z - k * dz * cos(theta); 
end


% This function calculates the directivity D for a linear broadside array antenna
function D = calculate_directivity_from_D_for_broadside(lambda, N, d)
    D = 2*N*d / (lambda);
end


% This function calculates the directivity from HPBW values
function D_HPBW = calculate_directivity_from_HPBW(theta_0_tonos, phi_0_tonos, Theta_x_tonos, Theta_y_tonos)
    Theta_h = 1 / (cos(theta_0_tonos) * sqrt((cos(phi_0_tonos)^2) / (Theta_x_tonos^2) + (sin(phi_0_tonos)^2) / (Theta_y_tonos^2)));
    Psi_h = 1 / (sqrt((cos(phi_0_tonos)^2) / (Theta_y_tonos^2) + (sin(phi_0_tonos)^2) / (Theta_x_tonos^2)));
    
    D_HPBW = 32400 / (Theta_h * Psi_h);
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

% Function to print the old and new theta and phi
function print_old_and_new_angles(theta_m, phi_m, theta_0_tonos, phi_0_tonos)
    fprintf('Old and New Angles (θ, φ) before and after the coordinate system swap:\n');
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
            
            % Print the result with deviations and broadside status
            fprintf('Old (θ, φ) = (%.2f, %.2f)° ---> ', ...
                    rad2deg(theta_old), rad2deg(phi_old));
            fprintf('New (θ, φ) = (%.2f, %.2f)°\n', ...
                    rad2deg(theta_new), rad2deg(phi_new));
        end
    end
end
