%% Antennas and Propagation Project - Nov 2025
% KARATIS DIMITRIOS 10775 
%% 1.4b

clc; close all; clear;   
%% Parameters and Constants (From 1.4a)
c0 = 3e8;                % Speed of light (m/s)
f = 1e9;                 % Frequency (1 GHz)
lambda = c0 / f;         % Wavelength
h0 = 120 * pi;           % Intrinsic impedance of free space (Ohms)
l = lambda / 2;          % Dipole length (l = lambda/2)
k = 2 * pi / lambda;     % Wave number
r = 50;                  % Distance to observation point in meters

%% Array Configuration
d_spacing = 3 * lambda / 4;  % Element spacing (d = 0.75 * lambda)
d_long = 2 * d_spacing;      % Distance between 1 and 3 (2d = 1.5 * lambda)

% Standard Self-Impedance Z11 (Half-Wave Dipole) -> Corresponds to 'a'
Z11_std = 73.1 + 1i * 42.5; 

%% Calculate Mutual Impedances Zmn (Z12, Z13) using Z_mutual (from 1.4a logic)
Z_mutual = @(d_input) calculate_mutual_impedance(d_input, k, l, h0);

% Z12 = Z(d=3λ/4) -> Corresponds to 'b'
Z12_mutual = Z_mutual(d_spacing); 

% Z13 = Z(d=3λ/2) -> Corresponds to 'c'
Z13_mutual = Z_mutual(d_long);

a = Z11_std;      % Z11 = a
b = Z12_mutual;   % Z12 = b
c = Z13_mutual;   % Z13 = c

%% Calculate Currents Ratios (I1/I2 = I3/I2)
% From the passive element equation (V1 = 0): a*I1 + b*I2 + c*I3 = 0
% Since I1 = I3: (a + c)*I1 = -b*I2
% Ratio I1/I2 (R_ratio) = -b / (a + c)
R_ratio = -b / (a + c); 

% Set current of the driven element (Element 2)
I2 = 1.0; 
I1 = R_ratio * I2; 
I3 = I1;

%% Calculate Input Impedance Z_in (Element 2)
% Z_in = V2 / I2
% V2 = b*I1 + a*I2 + b*I3
% Z_in = b*(I1/I2) + a + b*(I3/I2) = a + 2*b*(I1/I2)
Z_in = a + 2 * b * R_ratio; 

fprintf('1.4b\nCalculated System Parameters (using EMF method):\n');
fprintf('------------------------------------------------\n');
fprintf('  a = Z11 = %.2f + %.2fj Ohms\n', real(a), imag(a));
fprintf('  b = Z12 = %.2f + %.2fj Ohms\n', real(b), imag(b));
fprintf('  c = Z13 = %.2f + %.2fj Ohms\n', real(c), imag(c));
fprintf('\n');
fprintf('Current Ratio (I1/I2) = (I3/I2) = %.3f < %.1f deg\n', abs(R_ratio), rad2deg(angle(R_ratio)));
fprintf('Calculated Input Impedance Z_in = %.2f + %.2fj Ohms\n', real(Z_in), imag(Z_in));

%% Radiation Pattern Calculation (Horizontal Plane: theta = 90 deg)
theta_hor = pi/2; 
phi_hor = linspace(0, 2*pi, 360);

E0_abs_hor = abs(h0 * I2 * k * l * sin(theta_hor) / (4*pi*r));
        
% Array Factor
AF_hor = I2 + I1 .* exp(-1i * k * d_spacing * cos(phi_hor) * sin(theta_hor)) + I3 .* exp(1i * k * d_spacing * cos(phi_hor) * sin(theta_hor));
AF_abs_hor = abs(AF_hor);

% Normalize E
E_abs_hor = AF_abs_hor .* E0_abs_hor;
E_abs_hor_norm = E_abs_hor / max(E_abs_hor);

figure;
polarplot(phi_hor, E_abs_hor_norm, 'LineWidth', 2, 'Color', 'r');
title('Horizontal Radiation Pattern, d = 0.75\lambda');
ax = gca;
ax.ThetaDir = 'counterclockwise';
ax.ThetaZeroLocation = 'right'; 

 

%% Helper Function for Mutual Impedance Zmn (Logic from 1.4a)
function Z_mutual = calculate_mutual_impedance(d_input, k, l, h0)
    % Calculate the arguments U0, U1, U2
    U0 = k * d_input;
    U1 = k * (sqrt(d_input^2 + l^2) + l);
    U2 = k * (sqrt(d_input^2 + l^2) - l);

    % Calculate Sine and Cosine Integral Functions
    Si_U0 = sinint(U0);
    Si_U1 = sinint(U1);
    Si_U2 = sinint(U2);
    Ci_U0 = cosint(U0);
    Ci_U1 = cosint(U1);
    Ci_U2 = cosint(U2);

    % Calculate Mutual Impedance Zmn = Rmn + jXmn
    Rmn = (h0 / (4 * pi)) * (2 * Ci_U0 - Ci_U1 - Ci_U2);
    Xmn = -(h0 / (4 * pi)) * (2 * Si_U0 - Si_U1 - Si_U2);

    Z_mutual = Rmn + 1i * Xmn;
end