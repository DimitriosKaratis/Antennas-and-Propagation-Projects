%% Antennas and Propagation Project - Nov 2025
% KARATIS DIMITRIOS 10775 
%% 1.4a

clc; close all; clear;   
%% Parameters and Constants
c0 = 3e8;                % Speed of light (m/s)
f = 1e9;                 % Frequency (1 GHz)
lambda = c0 / f;         % Wavelength
h0 = 120 * pi;           % Intrinsic impedance of free space (Ohms)
l = lambda / 2;          % Dipole length (l = lambda/2)
k = 2 * pi / lambda;     % Wave number

%% Define Distance Range (d/lambda)
D_lambda = linspace(0.00001, 3.0, 300); 
d_meters = D_lambda * lambda;        % Convert to meters

%% Calculate Arguments U0, U1, U2
U0 = k .* d_meters;
U1 = k .* (sqrt(d_meters.^2 + l^2) + l);
U2 = k .* (sqrt(d_meters.^2 + l^2) - l);

%% Calculate Sine and Cosine Integral Functions Si(x) and Ci(x)
% Using MATLAB's built-in functions: sinint and cosint
Si_U0 = sinint(U0);
Si_U1 = sinint(U1);
Si_U2 = sinint(U2);

Ci_U0 = cosint(U0);
Ci_U1 = cosint(U1);
Ci_U2 = cosint(U2);

%% Calculate Mutual Impedance Z21m = R21m + jX21m
R21m = (h0 / (4 * pi)) .* (2 * Ci_U0 - Ci_U1 - Ci_U2);
X21m = -(h0 / (4 * pi)) .* (2 * Si_U0 - Si_U1 - Si_U2);

Z21m = R21m + 1i * X21m;

%% Plotting
figure;

% Plotting R21m (Real Part)
plot(D_lambda, R21m, 'LineWidth', 2, 'Color', [0 0.4470 0.7410], 'DisplayName', 'R_{21m} (Real Part)');
hold on;

% Plotting X21m (Imaginary Part)
plot(D_lambda, X21m, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'X_{21m} (Imaginary Part)');

% --- Adding Dashed Line at Z=0 ---
plot(D_lambda, zeros(size(D_lambda)), 'k--', 'LineWidth', 1, 'DisplayName', 'Z = 0 \Omega Line');

grid on;
title('Mutual Impedance Z_{21m} of Parallel \lambda/2 Dipoles (EMF Method)');
xlabel('Distance d / \lambda');
ylabel('Mutual Impedance Z_{21m} (Ohms)');
legend('show', 'Location', 'NorthEast');
xlim([0 3]);
ylim([-60 100]); 
