%% Antennas and Propagation Project - Nov 2025
% KARATIS DIMITRIOS 10775 
%% 1.3b

clc; clear; close all;

%% Parameters and Constants
N = 10;               % Number of elements
lambda = 1;           % Wavelength (in meters)
d = lambda / 2;       % Element spacing (half-wavelength)
k = 2 * pi / lambda;  % Wave number
delta = 0;            % Broadside configuration (delta = 0)

%% For SLL = -20 dBi
% Calculate directivity for non-uniform distribution
p_20 = [5.81259388666593; 5.47695152090878; 7.45631344887399; 8.24518675510177];
D_20 = calculate_D_full(p_20); 
D_20_dbi = 10 * log10(D_20); 
disp('1.3b');
disp('Directivity Calculations for SLL = -20 dB:');
disp('------------------------------------------');
disp(['Directivity D (Non-uniform, -20 dB) using formula from 1.2: ', num2str(D_20_dbi),' dBi']);

% Calculate directivity for Uniform distribution from the same formula
p_uni = [1.0, 1.0, 1.0, 1.0]; 
D_20_uni = calculate_D_full(p_uni); 
D_20_uni_dbi = 10 * log10(D_20_uni); 
disp(['Directivity D (Uniform, -20 dB) using formula from 1.2: ', num2str(D_20_uni_dbi),' dBi']);

% Calculate directivity using the formula from 1.1 (for uniform distribution)
D = calculate_directivity_from_D(k, N, d, delta);
D_dbi = 10 * log10(D); 
disp(['Directivity D (Uniform) using formula from 1.1: ', num2str(D_dbi),' dBi']);




%% For SLL = -30 dBi
disp(' '); disp(' ');
disp('Directivity Calculations for SLL = -30 dB:');
disp('------------------------------------------');

p_30 = [3.83493151254392; 5.63935063709595; 8.19052536606221; 9.47725214476735];
D_30 = calculate_D_full(p_30);
D_30_dbi = 10 * log10(D_30);
disp(['Directivity D (Non-uniform, -30 dB) using formula from 1.2: ', num2str(D_30_dbi),' dBi']);

% Calculate directivity for Uniform distribution from the same formula
p_uni = [1.0, 1.0, 1.0, 1.0]; 
D_30_uni = calculate_D_full(p_uni); 
D_30_uni_dbi = 10 * log10(D_30_uni); 
disp(['Directivity D (Uniform, -30 dB) using formula from 1.2: ', num2str(D_30_uni_dbi),' dBi']);

% Calculate directivity using the formula from 1.1 (for uniform distribution)
D = calculate_directivity_from_D(k, N, d, delta);
D_dbi = 10 * log10(D); 
disp(['Directivity D (Uniform) using formula from 1.1: ', num2str(D_dbi),' dBi']);




%% For SLL = -40 dBi
disp(' '); disp(' ');
disp('Directivity Calculations for SLL = -40 dB:');
disp('------------------------------------------');

p_40 = [2.56646139627147; 5.14168933000191; 8.01090353428315; 9.94538014324952];
D_40 = calculate_D_full(p_40); 
D_40_dbi = 10 * log10(D_40); 
disp(['Directivity D (Non-uniform, -40 dB) using formula from 1.2: ', num2str(D_40_dbi),' dBi']);

% Calculate directivity for Uniform distribution from the same formula
p_uni = [1.0, 1.0, 1.0, 1.0]; 
D_40_uni = calculate_D_full(p_uni); 
D_40_uni_dbi = 10 * log10(D_40_uni); 
disp(['Directivity D (Uniform, -40 dB) using formula from 1.2: ', num2str(D_40_uni_dbi),' dBi']);

% Calculate directivity using the formula from 1.1 (for uniform distribution)
D = calculate_directivity_from_D(k, N, d, delta);
D_dbi = 10 * log10(D); 
disp(['Directivity D (Uniform) using formula from 1.1: ', num2str(D_dbi),' dBi']);







%% HELPER FUNCTIONS

%  This function calculates the directivity D using the full formula
%  from Problem 1.2 for a non-uniform linear array antenna.
function D = calculate_D_full(p)

    % Parameter Setup
    N = 10;               % Number of elements
    lambda = 1.0;         % Wavelength (in meters)
    d = lambda / 2;       % Element spacing (half-wavelength)
    k = 2 * pi / lambda;  % Wave number
    delta = 0;           % Broadside configuration (delta = 0)
    
    % Current Distribution (Symmetric)
    I0 = 1.0; 
    % Symmetric current distribution: I = [I0, I1, I2, I3, I4, I4, I3, I2, I1, I0]
    I = [I0, p(1), p(2), p(3), p(4), p(4), p(3), p(2), p(1), I0];
    
    % Numerator Calculation
    Sum_I = sum(I);
    Numerator = k * d * (Sum_I)^2; 
    
    % Denominator Calculation
    Denominator = 0.0;
    
    % Loop over all elements to compute the denominator
    for n = 0:(N-1)
        for m = 0:(N-1)
            % Calculate the term I_n * I_m * exp(j(n-m)delta)
            I_term = I(n+1) * I(m+1) * exp(1j * (n - m) * delta);
            
            % Calculate sin[(n-m)kd] / (n-m)
            if n == m
                % Special case when n = m: The limit is kd
                Sin_kd_term = k * d;
            else
                % General case for n neq m
                n_minus_m_kd = (n - m) * k * d;
                n_minus_m = n - m;
                Sin_kd_term = sin(n_minus_m_kd) / n_minus_m;
            end
            
            % Add to the denominator
            Denominator = Denominator + I_term * Sin_kd_term;
        end
    end
    
    % Ensure the denominator is real (due to numerical errors)
    Denominator = real(Denominator);
    
    % Final Calculation formula
    D = Numerator / Denominator; 
end


% This function calculates the directivity D for a linear array antenna
function D = calculate_directivity_from_D(k, N, d , delta)
    % Define the function for Si(x), the integral of sin(t)/t from 0 to x
    Si = @(x) integral(@(t) sin(t) ./ t, 0, x);
    
    term1 = N * (-k * d + delta) / 2;  
    term2 = N * (k * d + delta) / 2;  
    
    if term2 < 10^(-5) 
        D = (N * k * d) / ((sin(term1)^2) / term1 - Si(term1 * 2));
    else
        D = (N * k * d) / ((sin(term1)^2) / term1 - (sin(term2)^2) / term2 +  Si(term2 * 2) - Si(term1 * 2));
    end
end