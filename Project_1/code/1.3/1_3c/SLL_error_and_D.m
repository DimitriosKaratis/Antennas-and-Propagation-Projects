%% Antennas and Propagation Project - Nov 2025
% KARATIS DIMITRIOS 10775 
%% Function for 1.3c

% SLL_error_and_D: Calculates two objectives for multi-objective optimization:
% 1. The mean squared error of the side lobes (SLL_error).
% 2. The negative directivity (-D) in linear scale.
% Both objectives are minimized.
function [f] = SLL_error_and_D(p, SLL_level_dB)
    
    N = 10;              % Number of elements
    lambda = 1.0;        % Wavelength
    d = lambda / 2;      % Element spacing
    delta = 0;           % Broadside configuration
    k = 2*pi / lambda;   % Wave number
    
    % Symmetric currents: I = [I0, I1, I2, I3, I4, I4, I3, I2, I1, I0]
    I0 = 1.0;
    I = [I0, p(1), p(2), p(3), p(4), p(4), p(3), p(2), p(1), I0]; 
    
    %% OBJECTIVE 1: SLL ERROR CALCULATION (MSE) 
    
    SLL_target = 10^(SLL_level_dB/20);
    
    % Define theta range and calculate psi
    theta = linspace(0, 90, 91);    
    psi = k*d*cosd(theta) + delta;  
    AF = zeros(1, length(theta)); 
    
    % Compute the Array Factor (AF)
    for n = 1:N
        AF = AF + I(n) * exp(1j*(n-1)*psi); 
    end
    
    % Normalized magnitude of the array factor |AF|
    AF_mag_norm = abs(AF) / max(abs(AF)); 
    [peaks, ~] = findpeaks(AF_mag_norm);
    
    % Remove the main lobe (absolute maximum)
    [max_peak, max_loc] = max(peaks);
    if max_peak > 0.98 
        peaks(max_loc) = []; 
    end
    
    if isempty(peaks)
        error_value = 1e6; % Large error if no side lobes are found
    else
        % (SLL_i - SLL_target)^2
        squared_error = (peaks - SLL_target).^2; 
        error_value = mean(squared_error); % Mean squared error (MSE)
    end
    

    %% OBJECTIVE 2: DIRECTIVITY CALCULATION (D)
    
    % Numerator Calculation (using formula from 1.2)
    Sum_I = sum(I);
    Numerator = k * d * (Sum_I)^2; 
    
    % Denominator Calculation
    Denominator = 0.0;
    
    for n = 0:(N-1)
        for m = 0:(N-1)
            % I_term = I_n * I_m * exp(j(n-m)delta)
            I_term = I(n+1) * I(m+1) * exp(1j * (n - m) * delta); % delta=0, so exp(0)=1
            
            if n == m
                Sin_kd_term = k * d; % Limit is k*d = pi
            else
                n_minus_m_kd = (n - m) * k * d;
                n_minus_m = n - m;
                Sin_kd_term = sin(n_minus_m_kd) / n_minus_m;
            end
            
            Denominator = Denominator + I_term * Sin_kd_term;
        end
    end
    
    % Ensure the denominator is real (due to numerical errors)
    D = Numerator / real(Denominator);
    
    %% RETURN OBJECTIVES (Minimize both)
    f(1) = error_value; % Minimize SLL_error
    f(2) = -D;          % Minimize -D (which maximizes D)
end