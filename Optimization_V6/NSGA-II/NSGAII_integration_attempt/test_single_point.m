% test_single_point.m - Script to test a single design point in isolation
% to diagnose GCI calculation issues
clear; clc; close all;

fprintf('==================================================\n');
fprintf('Testing single design point to diagnose GCI issues\n');
fprintf('==================================================\n\n');

% Define a test design point (conservative values)
base_rad = 0.12;           % Base platform radius (m)
top_rad = 0.05;            % Mobile platform radius (m)
h = 0.25;                  % Height (m)
theta_deg = 15;            % Pivot angle (degrees)
W_init = 0.01;             % Link width (m)
T_init = 0.01;             % Link thickness (m)
actuator_idx = 1;          % Actuator selection index
load_val = 9.81;           % Applied load (N·m)

fprintf('Test parameters:\n');
fprintf('  base_rad  = %.4f m\n', base_rad);
fprintf('  top_rad   = %.4f m\n', top_rad);
fprintf('  h         = %.4f m\n', h);
fprintf('  theta_deg = %.4f degrees\n', theta_deg);
fprintf('  W_init    = %.4f m\n', W_init);
fprintf('  T_init    = %.4f m\n', T_init);
fprintf('\n');

% Calculate alpha1 directly
alpha1 = atand(base_rad/h);
fprintf('Calculated alpha1 = %.4f degrees\n', alpha1);

% Leave alpha2 as empty to be calculated internally
alpha2 = [];

% Step 1: Run deflection optimization alone
fprintf('\n==== Step 1: Running deflection optimization ====\n');
debug_flag = true;  % Enable detailed output
[link_mass, isValid, W_opt, T_opt] = Deflection_Optim(load_val, h, theta_deg, 2*top_rad, debug_flag, W_init, T_init);
fprintf('Deflection optimization results:\n');
fprintf('  link_mass = %.2f g\n', link_mass);
fprintf('  isValid   = %d\n', isValid);
fprintf('  W_opt     = %.6f m\n', W_opt);
fprintf('  T_opt     = %.6f m\n', T_opt);

% Step 2: Calculate GCI directly
fprintf('\n==== Step 2: Calculating GCI directly ====\n');
GCI = SPM_kinematics_silent(base_rad, top_rad, h, theta_deg, alpha1, alpha2, isValid);
fprintf('Direct GCI calculation result: %.6f\n', GCI);

% Step 3: Test a manual grid of phi/theta/psi values
fprintf('\n==== Step 3: Testing sample workspace point directly ====\n');
phi_test = 0;
theta_test = 0;
psi_test = 0;
eta11 = 0; eta12 = 120; eta13 = 240;
beta2 = 90;
gamma2 = 120;

fprintf('Testing inverse kinematics for orientation:\n');
fprintf('  phi = %d, theta = %d, psi = %d\n', phi_test, theta_test, psi_test);

% Try to compute inverse kinematics for this orientation
try
    T = inverse_kinematics(alpha1, eta11, eta12, eta13, alpha2, beta2, phi_test, theta_test, psi_test);
    if any(isnan(T))
        fprintf('Inverse kinematics returned NaN values: [%.4f, %.4f, %.4f]\n', T(1), T(2), T(3));
    else
        fprintf('Joint angles from inverse kinematics: [%.4f, %.4f, %.4f] degrees\n', T(1), T(2), T(3));
        
        % Try forward kinematics
        [v1, v2, v3, w1, w2, w3] = forward_kinematics(alpha1, eta11, eta12, eta13, alpha2, gamma2, T(1), T(2), T(3));
        fprintf('Forward kinematics successful.\n');
        
        % Try calculating Jacobian
        J = calculate_Jacobian(w1, w2, w3, v1, v2, v3);
        fprintf('Jacobian calculation successful.\n');
        fprintf('Jacobian determinant: %.6e\n', det(J));
        
        % Try calculating conditioning index
        zeta = conditioning_index(J);
        fprintf('Conditioning index = %.6f\n', zeta);
    end
catch e
    fprintf('Error in calculation: %s\n', e.message);
end

% Step 4: Check multiple parameter combinations
fprintf('\n==== Step 4: Testing multiple parameter combinations ====\n');

% Define a range of parameter values to test
h_values = [0.2, 0.25, 0.3];
theta_values = [5, 10, 15];
base_rad_values = [0.08, 0.12, 0.16];
top_rad_values = [0.03, 0.05, 0.07];

% Initialize results matrix
results = zeros(length(h_values) * length(theta_values) * length(base_rad_values) * length(top_rad_values), 5);
result_idx = 1;

% Test all combinations
for h_idx = 1:length(h_values)
    for theta_idx = 1:length(theta_values)
        for base_idx = 1:length(base_rad_values)
            for top_idx = 1:length(top_rad_values)
                h_test = h_values(h_idx);
                theta_test = theta_values(theta_idx);
                base_test = base_rad_values(base_idx);
                top_test = top_rad_values(top_idx);
                
                % Calculate alpha1
                alpha1_test = atand(base_test/h_test);
                
                % Use empty alpha2 to be calculated internally
                alpha2_test = [];
                
                % Calculate GCI
                try
                    GCI_test = SPM_kinematics_silent(base_test, top_test, h_test, theta_test, alpha1_test, alpha2_test, true);
                    
                    % Store results
                    results(result_idx, :) = [h_test, theta_test, base_test, top_test, GCI_test];
                    result_idx = result_idx + 1;
                    
                    fprintf('h=%.2f, theta=%.1f, base=%.2f, top=%.2f: GCI=%.6f\n', h_test, theta_test, base_test, top_test, GCI_test);
                catch e
                    fprintf('Error for h=%.2f, theta=%.1f, base=%.2f, top=%.2f: %s\n', h_test, theta_test, base_test, top_test, e.message);
                end
            end
        end
    end
end

% Find best parameter combination
valid_results = results(results(:,5) > 0.001, :);
if ~isempty(valid_results)
    [max_gci, max_idx] = max(valid_results(:,5));
    best_params = valid_results(max_idx, 1:4);
    
    fprintf('\nBest parameter combination found:\n');
    fprintf('  h         = %.2f m\n', best_params(1));
    fprintf('  theta_deg = %.1f degrees\n', best_params(2));
    fprintf('  base_rad  = %.2f m\n', best_params(3));
    fprintf('  top_rad   = %.2f m\n', best_params(4));
    fprintf('  GCI       = %.6f\n', max_gci);
else
    fprintf('\nNo valid parameter combinations found with GCI > 0.001\n');
end

fprintf('\n==== Test complete ====\n');