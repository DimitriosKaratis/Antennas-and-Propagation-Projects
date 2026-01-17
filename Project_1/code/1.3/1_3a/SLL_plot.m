%% Antennas and Propagation Project - Nov 2025
% KARATIS DIMITRIOS 10775 
%% Function for 1.3a

% SLL_plot: Plots the normalized array factor in dB 
% for the optimal current distribution p.
function SLL_plot(p, SLL_target_dB)

    N = 10;              % Number of elements
    lambda = 1;          % Wavelength
    d = lambda / 2;      % d/lambda = 0.5 (half-wavelength spacing)
    delta = 0;           % delta = 0 (Broadside configuration)
    k = 2*pi / lambda;   % Wave number
    
    % Current distribution I = [I0, I1, I2, ..., I9, I0]
    I = [1.0, p(1), p(2), p(3), p(4), p(4), p(3), p(2), p(1), 1.0]; 

    % Compute the Array Factor (AF)
    theta = linspace(0, 90, 901);    
    psi = k*d*cosd(theta) + delta; 

    AF = zeros(1, length(theta)); 
    for n = 1:N
        AF = AF + I(n) * exp(1j*(n-1)*psi); 
    end
    
    % Normalization and Conversion to dB
    AF_mag_norm = abs(AF) / max(abs(AF)); % Normalize the magnitude of AF
    AF_dB = 20*log10(AF_mag_norm);        % Convert the normalized AF to dB
    
    % Plot the AF in dB
    plot(theta, AF_dB, 'LineWidth', 2); 
    hold on;
    
    % Plot the SLL target line
    line([0 90], [SLL_target_dB SLL_target_dB], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    
    % Identify side lobes (for visual confirmation)
    SLL_target_linear = 10^(SLL_target_dB/20); 
    [peaks, ~] = findpeaks(AF_mag_norm, 'MinPeakHeight', SLL_target_linear*0.00000001); 

    % Plot the side lobe levels
    for k = 1:length(peaks)
        % Find the angles corresponding to the k-th peak
        angle_peak = theta(AF_mag_norm == peaks(k));
        if isempty(angle_peak), continue; end
        
        peak_dB = 20*log10(peaks(k)); 
        
        % Mark only the side lobes (not the main lobe)
        if peak_dB < -1 % Exclude the main lobe (which is at 0 dB)
            scatter(angle_peak, peak_dB, 70, 'filled', 'MarkerFaceColor', 'b');
        end
    end
    
    grid on;
    xlabel('Angle \theta (degrees)');
    ylabel('Normalized Array Factor (dB)');
    title(sprintf('Radiation Pattern for SLL = %.1f dB', SLL_target_dB));
    legend('Array Factor', 'SLL Target');
    ylim([-60 0]); % Set lower limit for better visualization
    xlim([0 90]);
    hold off;
end
