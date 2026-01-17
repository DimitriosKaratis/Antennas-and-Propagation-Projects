%% Antennas and Propagation Project - Nov 2025
% KARATIS DIMITRIOS 10775 
%% 1.4c

clc; clear; close all;
%% Parameters and Constants (From 1.4a)
c0 = 3e8;                   % Speed of light (m/s)
f = 1e9;                    % Frequency (1 GHz)
lambda = c0 / f;            % Wavelength
h0 = 120 * pi;              % Intrinsic impedance of free space (Ohms)
l = lambda / 2;             % Dipole length (l = lambda/2)
k = 2 * pi / lambda;        % Wave number
Z11_std = 73.1 + 1i * 42.5; % Self-Impedance 'a'
Z0 = 50.0 + 0i;             % Characteristic Impedance Z0 = 50 Ohms

%% Define Distance Range (d/λ from 0 to 1)
% Range d in [0, λ]
D_lambda_range = linspace(0.001, 1.0, 400); 
d_meters_range = D_lambda_range * lambda;

%% Calculate Z_in(d)
% We use the helper function Z_mutual (from 1.4a logic) defined below.
Z_mutual = @(d_input) calculate_mutual_impedance(d_input, k, l, h0);

% b(d) = Z12(d)
b = Z_mutual(d_meters_range); 
% c(d) = Z13(2d)
c = Z_mutual(2 * d_meters_range);

a = Z11_std;
R_ratio = -b ./ (a + c); 
Z_in = a + 2 .* b .* R_ratio; 

%% Calculate Reflection Coefficient Gamma(d)
Gamma = (Z_in - Z0) ./ (Z_in + Z0);
Gamma_magnitude = abs(Gamma);

%% Plotting Magnitude of Gamma
figure;
plot(D_lambda_range, Gamma_magnitude, 'LineWidth', 2, 'Color', 'blue', 'DisplayName', '|Γ(d/λ)|');
hold on;
title('Reflection Coefficient Magnitude |Γ| vs. Distance d/λ (Parasitic Array, Z_0 = 50 Ohms)');
xlabel('Distance d / λ');
ylabel('Reflection Coefficient Magnitude |Γ|');
xlim([0 1.0]);
ylim([0 1.0]);
legend('Location', 'NorthEast');
grid on;

% Highlight Gamma = 0.3 threshold (Matching Criteria)
h_thresh = plot(D_lambda_range, 0.3 * ones(size(D_lambda_range)), 'r--', 'LineWidth', 1, 'DisplayName', '|Γ| = 0.3 Threshold');
uistack(h_thresh, 'bottom');

%% Analyze and print the matching regions (|Γ| < 0.3)
fprintf('1.4c\nReflection Coefficient Analysis (for d):\n');
fprintf('----------------------------------------\n');
[min_gamma, idx_min] = min(Gamma_magnitude);
fprintf('Minimum |Γ| = %.3f at d/λ = %.3f\n', min_gamma, D_lambda_range(idx_min));
fprintf('Regions where |Γ| < 0.3 (Good Matching):\n');

matching_region = Gamma_magnitude < 0.3;
indices = find(matching_region);
ranges = {};

% Logic to find contiguous segments where matching_region is True
if ~isempty(indices)
    current_start_index = indices(1);
    for i = 2:length(indices)
        % Check for a break in the sequence
        if indices(i) ~= indices(i-1) + 1
            ranges{end+1} = [D_lambda_range(current_start_index), D_lambda_range(indices(i-1))];
            current_start_index = indices(i);
        end
    end
    % Add the final segment
    ranges{end+1} = [D_lambda_range(current_start_index), D_lambda_range(indices(end))];
end

% Print the matching ranges
if isempty(ranges)
    fprintf('No distance d/λ in the range [0, 1] satisfies |Γ| < 0.3.\n');
else
    for range_pair = ranges
        start_d = range_pair{1}(1);
        end_d = range_pair{1}(2);
        fprintf('d/λ: %.3f to %.3f\n', start_d, end_d);
    end
end




%% Helper Function for Mutual Impedance Zmn (Logic from 1.4a)
function Z_mutual = calculate_mutual_impedance(d_input, k, l, h0)
    % Calculate the arguments U0, U1, U2
    U0 = k .* d_input;
    U1 = k .* (sqrt(d_input.^2 + l^2) + l);
    U2 = k .* (sqrt(d_input.^2 + l^2) - l);

    % Calculate Sine and Cosine Integral Functions
    Si_U0 = sinint(U0);
    Si_U1 = sinint(U1);
    Si_U2 = sinint(U2);
    Ci_U0 = cosint(U0);
    Ci_U1 = cosint(U1);
    Ci_U2 = cosint(U2);

    % Calculate Mutual Impedance Zmn = Rmn + jXmn
    Rmn = (h0 / (4 * pi)) .* (2 * Ci_U0 - Ci_U1 - Ci_U2);
    Xmn = -(h0 / (4 * pi)) .* (2 * Si_U0 - Si_U1 - Si_U2);

    Z_mutual = Rmn + 1i * Xmn;
end