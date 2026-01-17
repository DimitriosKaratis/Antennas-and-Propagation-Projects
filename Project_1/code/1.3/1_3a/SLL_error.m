%% Antennas and Propagation Project - Nov 2025
% KARATIS DIMITRIOS 10775 
%% Function for 1.3a

% SLL_error: Calculates the mean squared error of the side lobes 
% from the desired SLL (Side Lobe Level) in dB.
function [error_value] = SLL_error(p, SLL_level_dB)
  
    N = 10;              % Number of elements
    lambda = 1;          % Wavelength
    d = lambda / 2;      % d/lambda = 0.5 (half-wavelength spacing)
    delta = 0;           % delta = 0 (Broadside configuration)
    k = 2*pi / lambda;   % Wave number
    
    % Symmetric currents: I = [I0, I1, I2, ..., I9]
    % I0=I9=1.0, p = [I1, I2, I3, I4]
    I = [1.0, p(1), p(2), p(3), p(4), p(4), p(3), p(2), p(1), 1.0]; 

    % Desired SLL level in linear scale (from dB)
    SLL_target = 10^(SLL_level_dB/20);

    % Compute the Array Factor (AF)
    theta = linspace(0, 90, 91);    
    psi = k*d*cosd(theta) + delta;  

    AF = zeros(1, length(theta)); 
    for n = 1:N
        AF = AF + I(n) * exp(1j*(n-1)*psi); 
    end
    
    % Normalized magnitude of the array factor |AF|
    AF_mag_norm = abs(AF) / max(abs(AF)); 

    % Identify Side Lobes
    [peaks, ~] = findpeaks(AF_mag_norm);
    
    % Exclude the Main Lobe. The main lobe is the absolute maximum (1.0).
    % For Broadside (delta=0), the main lobe is at theta=90 (loc=91).
    [max_peak, max_loc] = max(peaks);
    
    % Remove the main lobe if it's within the range [0, 90]
    if max_peak > 0.98 
        peaks(max_loc) = []; 
    end

    % Compute the Error (Mean Squared Error)
    if isempty(peaks)
        % If no side lobes are found, return a large error value
        error_value = 1e6;
    else
        % (SLL_i - SLL_target)^2
        squared_error = (peaks - SLL_target).^2; 
        error_value = mean(squared_error); % Mean squared error (MSE)
    end
end
