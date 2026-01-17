%% Antennas and Propagation Project - Nov 2025
% KARATIS DIMITRIOS 10775 
%% 1.3a

clc; clear; close all;

% Create a new figure for the subplots
figure;

% First subplot: SLL Plot for -20 dB
SLL_target_dB = -20; 
p_20 = [5.81259388666593; 5.47695152090878; 7.45631344887399; 8.24518675510177];

subplot(3, 2, 1); 
SLL_plot(p_20, SLL_target_dB);
hold on;

% Second subplot: SLL Plot (AF vs. psi) for -20 dB
subplot(3, 2, 2); 
SLL_plot_psi(p_20, SLL_target_dB);
hold on;

% Third subplot: SLL Plot for -30 dB
SLL_target_dB = -30; 
p_30 = [3.83493151254392; 5.63935063709595; 8.19052536606221; 9.47725214476735];

subplot(3, 2, 3); 
SLL_plot(p_30, SLL_target_dB);
hold on;

% Fourth subplot: SLL Plot (AF vs. psi) for -30 dB
subplot(3, 2, 4); 
SLL_plot_psi(p_30, SLL_target_dB);
hold on;

% Fifth subplot: SLL Plot for -40 dB
SLL_target_dB = -40; 
p_40 = [2.56646139627147; 5.14168933000191; 8.01090353428315; 9.94538014324952];

subplot(3, 2, 5);
SLL_plot(p_40, SLL_target_dB);
hold on;

% Sixth subplot: SLL Plot (AF vs. psi) for -40 dB
subplot(3, 2, 6); 
SLL_plot_psi(p_40, SLL_target_dB);
hold off;



%% --- PLOTS FOR TESTING AND UNDERSTANDING SIDE-LOBE SUPPRESSION --- 
figure;

% First subplot: Plot for p = [1 1 1 1 1 1 1 1 1 1]
subplot(3, 2, 1); 
p = [1 1 1 1];
SLL_plot_psi_for_I(p);
hold on;

% Second subplot: Plot for p = [1 1 1 1 2 2 1 1 1 1]
subplot(3, 2, 3); 
p = [1 1 1 2];
SLL_plot_psi_for_I(p);
hold on;

% Third subplot: Plot for p = [1 1 1 2 2 2 1 1 1 1]
subplot(3, 2, 5); 
p = [1 1 2 2];
SLL_plot_psi_for_I(p);
hold on;

% Fourth subplot: Plot for p = [1 1 2 2 2 2 2 2 1 1]
subplot(3, 2, 2); 
p = [1 2 2 2];
SLL_plot_psi_for_I(p);
hold on;

% Fifth subplot: Plot for p = [1 2 2 2 2 2 2 2 2 1]
subplot(3, 2, 4); 
p = [2 2 2 2];
SLL_plot_psi_for_I(p);
hold on;

% Sixth subplot: Plot for p = [1 2 2 2 3 3 2 2 2 1]
subplot(3, 2, 6); 
p = [2 2 2 3];
SLL_plot_psi_for_I(p);
hold off;






% Helper function that plots |AF| vs psi
function SLL_plot_psi(p, SLL_target_dB)
    N = 10;              % Number of elements
    lambda = 1;          % Wavelength
    d = lambda / 2;      % d/lambda = 0.5 (half-wavelength spacing)
    delta = 0;           % delta = 0 (Broadside configuration)
    k = 2*pi / lambda;   % Wave number
    
    % Current distribution I = [I0, I1, I2, ..., I9, I0]
    I = [1.0, p(1), p(2), p(3), p(4), p(4), p(3), p(2), p(1), 1.0]; 
    
    % Compute the Array Factor (AF)
    theta = linspace(-180, -90, 901); 
    psi = k*d*cosd(theta) + delta; 

    AF = zeros(1, length(theta)); 
    for n = 1:N
        AF = AF + I(n) * exp(1j*(n-1)*psi); 
    end
    
    % Normalization and Conversion to dB
    AF_mag_norm = abs(AF) / max(abs(AF)); 
        
    plot(psi, AF_mag_norm, 'LineWidth', 2); 
    title(sprintf('Radiation Pattern (|AF| vs. \\psi) for SLL = %.1f dB', SLL_target_dB));
    xlabel('Spatial Variable \psi (radians)');
    ylabel('Normalized Array Factor (dB)');
    grid on;
end


% Helper function that plots |AF| vs psi for given p values. 
% Used for testing and understanding patterns.
function SLL_plot_psi_for_I(p)
    N = 10;              % Number of elements
    lambda = 1;          % Wavelength
    d = lambda / 2;      % d/lambda = 0.5 (half-wavelength spacing)
    delta = 0;           % delta = 0 (Broadside configuration)
    k = 2*pi / lambda;   % Wave number
    
    % Current distribution I = [I0, I1, I2, ..., I9, I0]
    I = [1.0, p(1), p(2), p(3), p(4), p(4), p(3), p(2), p(1), 1.0]; 
    
    % Compute the Array Factor (AF)
    theta = linspace(-180, -90, 901); 
    psi = k*d*cosd(theta) + delta; 

    AF = zeros(1, length(theta)); 
    for n = 1:N
        AF = AF + I(n) * exp(1j*(n-1)*psi); 
    end
    
    % Normalization and Conversion to dB
    AF_mag_norm = abs(AF) / max(abs(AF)); 
        
    plot(psi, AF_mag_norm, 'LineWidth', 2); 
    title(sprintf(['Radiation Pattern (|AF| vs. \\psi) for' ...
        '\nI = [%.1d, %.1d, %.1d, %.1d, %.1d, %.1d, %.1d, %.1d, %.1d, %.1d]'], I));
    xlabel('Spatial Variable \psi (radians)');
    ylabel('Normalized Array Factor (dB)');
    grid on;
end



