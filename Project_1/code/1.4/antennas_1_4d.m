%% Antennas and Propagation Project - Nov 2025
% KARATIS DIMITRIOS 10775 
%% 1.4d

clc; clear; close all;
%% Parameters and Constants
c0 = 3e8;                % Speed of light (m/s)
f = 1e9;                 % Frequency (1 GHz)
lambda = c0 / f;         % Wavelength
h0 = 120 * pi;           % Intrinsic impedance of free space (Ohms)
l = lambda / 2;          % Dipole length (l = lambda/2)
k = 2 * pi / lambda;     % Wave number
Z0 = 50.0 + 0i;          % Characteristic Impedance Z0 = 50 Ohms

%% Define 2D Parameter Space
% Range d in [0, λ] (Spacing between elements)
D_lambda_range = linspace(0.01, 1.0, 50); 
% Range h in [0, λ] (Height above ground plane)
H_lambda_range = linspace(0.01, 1.0, 50); 

[D_grid, H_grid] = meshgrid(D_lambda_range, H_lambda_range);

d_meters_grid = D_grid * lambda;
h_meters_grid = H_grid * lambda;

%% Calculate Zmn Componentσ
% Standard Self-Impedance Z11 (Half-Wave Dipole, thin wire)
Z11 = 73.1 + 1i * 42.5; 
Z22 = Z11;

% Helper function Z_mutual (from 1.4a logic) defined below.
Z_mutual_parallel = @(d_input) calculate_mutual_impedance_parallel(d_input, k, l, h0);

Z12 = Z_mutual_parallel(d_meters_grid); 
Z21 = Z12;

Z13 = Z_mutual_parallel(2 * d_meters_grid);
Z23 = Z12;

Z14 = Z_mutual_parallel(sqrt((2*d_meters_grid).^2 + (2*h_meters_grid).^2));

Z15 = Z_mutual_parallel(sqrt(d_meters_grid.^2 + (2*h_meters_grid).^2));
Z24 = Z15;
Z26 = Z15;

Z16 = Z_mutual_parallel(2*h_meters_grid);
Z25 = Z16;

R_ratio = (Z15-Z12) ./ (Z11+Z13-Z14-Z16); 
Z_in = Z21.*R_ratio + Z22 + Z23.*R_ratio -Z24.*R_ratio -Z25 -Z26.*R_ratio; 




%% Calculate Reflection Coefficient Gamma(d, h)
Gamma_new = (Z_in- Z0) ./ (Z_in + Z0);
Gamma_magnitude = abs(Gamma_new);


%% Plotting 2D Contour 
figure;
contourf(D_grid, H_grid, Gamma_magnitude, [0 0.1 0.2 0.3 0.5 0.7 1.0], 'LineWidth', 0.5);
colorbar; 

hold on;
contour(D_grid, H_grid, Gamma_magnitude, [0.3 0.3], 'k-', 'LineWidth', 2);

title('Contour: |Γ| vs. d/\lambda and h/\lambda');
xlabel('Element Spacing d / \lambda');
ylabel('Height above Reflector h / \lambda');
xlim([0 1.0]);
ylim([0 1.0]);
colormap jet; 

%% Plotting 3D Surface
figure;
surf(D_grid, H_grid, Gamma_magnitude, Gamma_magnitude, 'EdgeColor', 'none');
view(3);

zlabel('Reflection Coeff |Γ|');
xlabel('Element Spacing d / \lambda');
ylabel('Height above Reflector h / \lambda');
title('Surface: |Γ| 3D Visualization');
colorbar;
zlim([0 1.0]); 

%% Find points where |Gamma| < 0.3
[rows, cols] = find(Gamma_magnitude < 0.3);  

matching_points = [];

% Loop through the found points and collect the corresponding values of d/λ and h/λ
for i = 1:length(rows)
    d_lambda = D_grid(rows(i), cols(i));  % d/λ at the matching point
    h_lambda = H_grid(rows(i), cols(i));  % h/λ at the matching point
    matching_points = [matching_points; d_lambda, h_lambda]; 
end

% Print the results
fprintf('1.4d\nReflection Coefficient Analysis (for d and h): \n');
fprintf('----------------------------------------------\n');
min_gamma = min(Gamma_magnitude(:));
[row, col] = find(Gamma_magnitude == min_gamma, 1);
fprintf('Minimum |Γ| = %.3f at d/λ = %.3f, h/λ = %.3f\n', min_gamma, D_grid(row, col), H_grid(row, col));

fprintf('\nMatching Points where |Γ| < 0.3:\n');
fprintf('\n');
fprintf('d/λ       h/λ\n');
fprintf('---------- --------\n');
for i = 1:size(matching_points, 1)
    fprintf('%.3f      %.3f\n', matching_points(i, 1), matching_points(i, 2));
end







%% Helper Function for Mutual Impedance Zmn (Logic from 1.4a)
function Z_mutual_parallel = calculate_mutual_impedance_parallel(d_input, k, l, h0)
    % Calculate the arguments U0, U1, U2
    U0 = k .* d_input;
    U1 = k .* (sqrt(d_input.^2 + l^2) + l);
    U2 = k .* (sqrt(d_input.^2 + l^2) - l);

    % Calculate Sine and Cosine Integral Functions
    Si_U0 = sinint(U0);
    Ci_U0 = cosint(U0);
    Si_U1 = sinint(U1);
    Ci_U1 = cosint(U1);
    Si_U2 = sinint(U2);
    Ci_U2 = cosint(U2);

    % Calculate Mutual Impedance Zmn = Rmn + jXmn
    Rmn = (h0 / (4 * pi)) .* (2 * Ci_U0 - Ci_U1 - Ci_U2);
    Xmn = -(h0 / (4 * pi)) .* (2 * Si_U0 - Si_U1 - Si_U2);

    Z_mutual_parallel = Rmn + 1i * Xmn;
end
