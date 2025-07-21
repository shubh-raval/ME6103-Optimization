% quick_test_spm.m - Script to test optimization with the best parameters
clear; clc; close all;

fprintf('==================================================\n');
fprintf('Testing SPM optimization with best parameters\n');
fprintf('==================================================\n\n');

% Test the optimal parameter combination found earlier
params.load = 9.81;           % Applied load (N·m)
params.platform_diam = 0.10;  % Platform diameter (m) 
params.base_rad = 0.16;       % Base platform radius (m) - optimal from testing
params.h = 0.20;              % Height (m) - optimal from testing
params.theta_deg = 15;        % Pivot axis angle (deg) - optimal from testing
params.actuator_masses = [282, 107, 230, 260, 570, 195, 248, 248, 41, 90, 450, 450, 207.8, 120, 135, 300, 4000, 122];
params.debug = true;          % Enable detailed output

% Test with different link dimensions
W_values = [0.01, 0.02, 0.03];
T_values = [0.01, 0.02, 0.03];
actuator_idx = 9;  % Index 9 corresponds to the lightest actuator (41g)

% Create results array
results = zeros(length(W_values) * length(T_values), 4);
result_idx = 1;

fprintf('Testing with best geometric parameters:\n');
fprintf('  h = %.2f m\n', params.h);
fprintf('  theta_deg = %.1f degrees\n', params.theta_deg);
fprintf('  base_rad = %.2f m\n', params.base_rad);
fprintf('  actuator_idx = %d (mass = %.1f g)\n\n', actuator_idx, params.actuator_masses(actuator_idx));

% Test all W,T combinations
for i = 1:length(W_values)
    for j = 1:length(T_values)
        W_init = W_values(i);
        T_init = T_values(j);
        
        fprintf('Testing W=%.3f, T=%.3f:\n', W_init, T_init);
        
        % Evaluate objectives
        F = SPM_objectives([W_init, T_init, actuator_idx], params);
        
        % Store results [W, T, mass, GCI]
        results(result_idx, :) = [W_init, T_init, F(1), -F(2)];
        result_idx = result_idx + 1;
        
        fprintf('  Mass = %.2f g, GCI = %.6f\n\n', F(1), -F(2));
    end
end

% Find best solution
[~, best_idx] = max(results(:,4));
fprintf('Best solution from test:\n');
fprintf('  W = %.4f m\n', results(best_idx, 1));
fprintf('  T = %.4f m\n', results(best_idx, 2));
fprintf('  Mass = %.2f g\n', results(best_idx, 3));
fprintf('  GCI = %.6f\n', results(best_idx, 4));

% Visualize results
figure;
scatter3(results(:,1), results(:,2), results(:,4), 100, results(:,3), 'filled');
xlabel('Width (m)');
ylabel('Thickness (m)');
zlabel('GCI');
c = colorbar;
ylabel(c, 'Mass (g)');
title('SPM Design Space (Fixed h=0.2m, \theta=15°)');
grid on;

fprintf('\nRecommended starting point for optimization:\n');
fprintf('  W = 0.02 m (or less for lighter mass)\n');
fprintf('  T = 0.02 m (or less for lighter mass)\n');
fprintf('  actuator_idx = 9 (41g actuator)\n');
fprintf('  h = 0.20 m\n');
fprintf('  theta_deg = 15 degrees\n');