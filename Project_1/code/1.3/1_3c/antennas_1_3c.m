%% Antennas and Propagation Project - Nov 2025
% KARATIS DIMITRIOS 10775 
% 1.3c

clc; clear; close all; 

% Define SLL values for the loop
SLL_values = [-20, -30, -40];


disp('1.3c');
% Loop through each SLL value
for i = 1:length(SLL_values)
    
    % Load the Pareto front results (SLL error and -D)
    FVAL_struct = load(['Pareto_objValue_SLL', num2str(SLL_values(i)*(-1)), '.mat']);
    Currents_struct = load(['Pareto_sol_SLL', num2str(SLL_values(i)*(-1)), '.mat']); 
    
    % Currents matrix (42 x 4): Current solutions
    Currents = Currents_struct.solution; 

    % FVAL matrix (42 x 2): Objective values
    FVAL = FVAL_struct.objectiveValue; 

    % Find Best Solution (Min SLL_error)
    % Find index (idx) for the row with the minimum SLL error (first column)
    [min_error_value, idx] = min(FVAL(:, 1)); 

    % Best currents [I1, I2, I3, I4]
    Optimal_Currents = Currents(idx, :); 

    % Calculate Directivity D (D = -Objective 2)
    Directivity_D = -FVAL(idx, 2); 

    % Display Results 
    disp('-----------------------------------------------------------');
    disp(['Best Pareto Solution for SLL = ', num2str(SLL_values(i)), ' dB (Min SLL Error):']);
    disp('-----------------------------------------------------------');
    disp(['1. Minimum SLL Error (Objective 1): ', num2str(min_error_value, '%0.2e')]);
    disp(['2. Maximum Directivity D (dBi):     ', num2str(10 * log10(Directivity_D), '%0.4f'), ' dBi']);
    disp('3. Optimal Currents [I1, I2, I3, I4]: ');
    disp(Optimal_Currents);

end
