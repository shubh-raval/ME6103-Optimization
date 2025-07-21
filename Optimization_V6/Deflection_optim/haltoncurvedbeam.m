% Load and process beam experiment data
clear all;
close all;

% Load the existing beam configurations
beam_config_table = readtable('beam_configurations.csv');
beam_params = [beam_config_table.a_m, beam_config_table.b_m, beam_config_table.w_m, beam_config_table.t_m];

% Define loads
loads = [100, 300, 500, 700]; % N

% Initialize results matrix [beam_number, a, b, w, t, load, deflection]
num_tests = length(loads) * 20;
results = zeros(num_tests, 7);

% Fill in known parameters
test_idx = 1;
for beam_idx = 1:20
    for load_idx = 1:length(loads)
        results(test_idx, 1) = beam_idx;
        results(test_idx, 2:5) = beam_params(beam_idx, :);
        results(test_idx, 6) = loads(load_idx);
        test_idx = test_idx + 1;
    end
end

% Input measured deflection values
deflections = zeros(80, 1); % Initialize with zeros

%% BEAM 1
deflections(1) = 6.123e-02;  % mm
deflections(2) = 1.839e-01;  % mm
deflections(3) = 3.070e-01;  % mm
deflections(4) = 4.304e-01;  % mm

%% BEAM 2
deflections(5) = 1.597e-01;  % mm
deflections(6) = 4.792e-01;  % mm
deflections(7) = 7.988e-01;  % mm
deflections(8) = 1.118;      % mm

%% BEAM 3
deflections(9)  = 2.118e-01; % mm
deflections(10) = 6.356e-01; % mm
deflections(11) = 1.06;      % mm
deflections(12) = 1.484;     % mm

%% BEAM 4
deflections(13) = 3.372e-03; % mm
deflections(14) = 1.012e-02; % mm
deflections(15) = 1.686e-02; % mm
deflections(16) = 2.361e-02; % mm

%% BEAM 5
deflections(17) = 9.253e-02; % mm
deflections(18) = 2.777e-01; % mm
deflections(19) = 4.631e-01; % mm
deflections(20) = 6.486e-01; % mm

%% BEAM 6
deflections(21) = 1.786e-01; % mm
deflections(22) = 5.36e-01;  % mm
deflections(23) = 8.936e-01; % mm
deflections(24) = 1.251;     % mm

%% BEAM 7
deflections(25) = 4.263e-02; % mm
deflections(26) = 1.279e-01; % mm
deflections(27) = 2.132e-01; % mm
deflections(28) = 2.985e-01; % mm

%% BEAM 8
deflections(29) = 3.258e-02; % mm
deflections(30) = 9.775e-02; % mm
deflections(31) = 1.629e-01; % mm
deflections(32) = 2.281e-01; % mm

%% BEAM 9
deflections(33) = 2.671;     % mm
deflections(34) = 8.022;     % mm
deflections(35) = 13.35;     % mm - Fixed from 1.335
deflections(36) = 18.61;     % mm - Fixed from 1.861e01

%% BEAM 10
deflections(37) = 2.678e-02; % mm
deflections(38) = 8.034e-02; % mm - Fixed from .8034e-02
deflections(39) = 1.339e-01; % mm
deflections(40) = 1.875e-01; % mm

%% BEAM 11
deflections(41) = 3.866e-02; % mm
deflections(42) = 1.16e-01;  % mm
deflections(43) = 1.934e-01; % mm
deflections(44) = 2.708e-01; % mm

%% BEAM 12
deflections(45) = 1.85e-01;  % mm
deflections(46) = 5.559e-01; % mm
deflections(47) = 9.266e-01; % mm
deflections(48) = 1.297;     % mm

%% BEAM 13
deflections(49) = 3.370e-02; % mm
deflections(50) = 1.012e-01; % mm
deflections(51) = 1.687e-01; % mm
deflections(52) = 2.364e-01; % mm

%% BEAM 14
deflections(53) = 2.354e-02; % mm
deflections(54) = 7.067e-02; % mm
deflections(55) = 1.179e-01; % mm
deflections(56) = 1.651e-01; % mm

%% BEAM 15
deflections(57) = 4.664e-01; % mm
deflections(58) = 1.4;       % mm
deflections(59) = 2.334;     % mm
deflections(60) = 3.268;     % mm

%% BEAM 16
deflections(61) = 2.442e-03; % mm
deflections(62) = 7.326e-03; % mm
deflections(63) = 1.221e-02; % mm
deflections(64) = 1.709e-02; % mm

%% BEAM 17
deflections(65) = 2.868e-01; % mm
deflections(66) = 8.621e-01; % mm
deflections(67) = 1.440;     % mm
deflections(68) = 2.019;     % mm

%% BEAM 18
deflections(69) = 1.150;     % mm
deflections(70) = 3.454;     % mm
deflections(71) = 5.758;     % mm
deflections(72) = 8.06;      % mm

%% BEAM 19
deflections(73) = 1.501e-02; % mm
deflections(74) = 4.503e-02; % mm
deflections(75) = 7.506e-02; % mm
deflections(76) = 1.501e-01; % mm

%% BEAM 20
deflections(77) = 6.38e-02;  % mm
deflections(78) = 1.916e-01; % mm
deflections(79) = 3.194e-01; % mm
deflections(80) = 4.472e-01; % mm

% Add deflection values to the results table
results(:,7) = deflections;
results_table = array2table(results, 'VariableNames', {'BeamNumber', 'a_m', 'b_m', 'w_m', 't_m', 'Load_N', 'Deflection_mm'});

% Save results
writetable(results_table, 'beam_test_results.csv');
disp('Experimental results saved to "beam_test_results.csv"');

%% Data visualization before model fitting
figure;
subplot(2,2,1);
scatter3(results(:,3), results(:,6), results(:,7), 30, results(:,1), 'filled');
xlabel('Width (m)');
ylabel('Load (N)');
zlabel('Deflection (mm)');
title('Deflection vs Width and Load');
cbar = colorbar;
cbar.Label.String = 'Beam #';

subplot(2,2,2);
scatter3(results(:,4), results(:,6), results(:,7), 30, results(:,1), 'filled');
xlabel('Thickness (m)');
ylabel('Load (N)');
zlabel('Deflection (mm)');
title('Deflection vs Thickness and Load');
cbar = colorbar;
cbar.Label.String = 'Beam #';

subplot(2,2,3);
scatter3(results(:,2), results(:,6), results(:,7), 30, results(:,1), 'filled');
xlabel('a (m)');
ylabel('Load (N)');
zlabel('Deflection (mm)');
title('Deflection vs a and Load');
cbar = colorbar;
cbar.Label.String = 'Beam #';

subplot(2,2,4);
boxplot(results(:,7), results(:,6));
xlabel('Load (N)');
ylabel('Deflection (mm)');
title('Deflection Distribution by Load');

%% Build meta-model (using Gaussian Process Regression)
X = results(:, 2:6); % Features: a, b, w, t, load
y = results(:, 7);   % Target: deflection

% Create a GPR model
disp('Training the Gaussian Process Regression model...');
gprMdl = fitrgp(X, y, 'KernelFunction', 'ardsquaredexponential', ...
                      'Standardize', true);

% Evaluate model on training data
ypred = predict(gprMdl, X);
rmse = sqrt(mean((y - ypred).^2));
rsquared = 1 - sum((y - ypred).^2)/sum((y - mean(y)).^2);

fprintf('Meta-model performance: RMSE = %.4f mm, R² = %.4f\n', rmse, rsquared);

% Plot actual vs predicted values
figure;
scatter(y, ypred, 40, results(:,1), 'filled');
hold on;
plot([min(y), max(y)], [min(y), max(y)], 'r--', 'LineWidth', 2);
xlabel('Actual Deflection (mm)');
ylabel('Predicted Deflection (mm)');
title('Meta-model Performance: Actual vs. Predicted');
cbar = colorbar;
cbar.Label.String = 'Beam #';
grid on;

% Save the model
save('beam_deflection_metamodel.mat', 'gprMdl');
disp('Meta-model saved to "beam_deflection_metamodel.mat"');

%% Model inference example
% Example prediction for a new beam
new_beam = [0.22, 0.07, 0.025, 0.015, 100]; % [a, b, w, t, load]
predicted_deflection = predict(gprMdl, new_beam);

% Cantilever Beam
a = new_beam(1);
b = new_beam(2);
P = new_beam(5);
w = new_beam(3);
t = new_beam(4);
I = (w * t^3) / 12; % second moment of area for rectangular cross-section
E = 69e9;          % Young's modulus in Pascals
quarter_circumference = (pi/4) * (3*(a + b) - sqrt((3*a + b)*(a + 3*b)));
deflection_cb = P * quarter_circumference^3 / (3 * E * I);

fprintf('Example prediction - Beam with a=%.4f, b=%.4f, w=%.4f, t=%.4f at %.0fN load\n', ...
        new_beam(1), new_beam(2), new_beam(3), new_beam(4), new_beam(5));
fprintf('Gaussian Metamodel deflection: %.4f mm\n', predicted_deflection);
fprintf('Cantilever Beam Equivalent deflection: %.4f mm\n', deflection_cb);
% Add a function to explore the parameter space
fprintf('\nYou can now use the trained model to predict deflection for any beam within your parameter space.\n');
fprintf('Example usage: predictDeflection(gprMdl, [a, b, w, t, load])\n');

% Parameter space exploration (add a figure to show sensitivity)
figure;
% Create a grid of loads
grid_loads = linspace(100, 700, 50);

% Pick a few representative beams
beam_indices = [4, 16, 10, 20, 9, 18]; % From low to high deflection

% Predict deflection across load range for these beams
for i = 1:length(beam_indices)
    beam_idx = beam_indices(i);
    beam_params_fixed = repmat(beam_params(beam_idx,:), length(grid_loads), 1);
    grid_X = [beam_params_fixed(:,1), beam_params_fixed(:,2), beam_params_fixed(:,3), beam_params_fixed(:,4), grid_loads'];
    grid_pred = predict(gprMdl, grid_X);
    
    % Plot load vs deflection curve
    plot(grid_loads, grid_pred, 'LineWidth', 2);
    hold on;
    
    % Add actual experimental points
    beam_results_idx = find(results(:,1) == beam_idx);
    scatter(results(beam_results_idx, 6), results(beam_results_idx, 7), 80, 'filled');
end

xlabel('Load (N)');
ylabel('Deflection (mm)');
title('Load vs Deflection for Different Beam Configurations');
legend({'Beam 4', 'Actual', 'Beam 16', 'Actual', 'Beam 10', 'Actual', ...
        'Beam 20', 'Actual', 'Beam 9', 'Actual', 'Beam 18', 'Actual'}, ...
        'Location', 'northwest');
grid on;

% Function for easy prediction
function deflection = predictDeflection(model, params)
    % params should be [a, b, w, t, load]
    deflection = predict(model, params);
    fprintf('Predicted deflection: %.4f mm\n', deflection);
end

%% Visualize the model vs the predicted across width and thickness
figure;
scatter3(results(:,3), results(:,6), results(:,7), 30, results(:,1), 'filled');
xlabel('Width (m)');
ylabel('Load (N)');
zlabel('Deflection (mm)');
title('Deflection vs Width and Load');
cbar = colorbar;
cbar.Label.String = 'Beam #';
