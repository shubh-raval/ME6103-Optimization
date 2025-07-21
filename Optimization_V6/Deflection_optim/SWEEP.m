% Parameter Sweep Script for Deflection_Optim
% This script evaluates the optimization function across a range of inputs
% to understand sensitivity and performance characteristics

% Create storage arrays for results
results = struct();

% 1. LOAD SWEEP - Testing different torque values
fprintf('\n===== LOAD SWEEP =====\n');
load_range = [2.5, 5, 7.5, 10, 15, 20];  % Range of torque values (N·m)
results.load_sweep.loads = load_range;
results.load_sweep.mass = zeros(size(load_range));
results.load_sweep.valid = zeros(size(load_range));

% Fixed parameters for load sweep
H_val = 0.15;      % Mechanism height (m)
theta_deg = 5;     % Pivot axis angle (degrees)
platform_diam = 0.14;  % Platform diameter (m)

for i = 1:length(load_range)
    fprintf('\nTesting Load = %.2f N·m\n', load_range(i));
    % Run with debug=false to suppress detailed output for each run
    [mass, isValid] = Deflection_Optim(load_range(i), H_val, theta_deg, platform_diam, false);
    results.load_sweep.mass(i) = mass;
    results.load_sweep.valid(i) = isValid;
    fprintf('  Mass: %.2f g, Valid: %d\n', mass, isValid);
end

% 2. ANGLE SWEEP - Testing different pivot angles
fprintf('\n===== ANGLE SWEEP =====\n');
angle_range = [0, 5, 10, 15, 20, 30, 45];  % Range of pivot angles (degrees)
results.angle_sweep.angles = angle_range;
results.angle_sweep.mass = zeros(size(angle_range));
results.angle_sweep.valid = zeros(size(angle_range));

% Fixed parameters for angle sweep
load_val = 9.81;   % Applied load (N·m)
H_val = 0.15;      % Mechanism height (m)
platform_diam = 0.14;  % Platform diameter (m)

for i = 1:length(angle_range)
    fprintf('\nTesting Angle = %.2f degrees\n', angle_range(i));
    [mass, isValid] = Deflection_Optim(load_val, H_val, angle_range(i), platform_diam, false);
    results.angle_sweep.mass(i) = mass;
    results.angle_sweep.valid(i) = isValid;
    fprintf('  Mass: %.2f g, Valid: %d\n', mass, isValid);
end

% 3. HEIGHT SWEEP - Testing different mechanism heights
fprintf('\n===== HEIGHT SWEEP =====\n');
height_range = [0.1, 0.125, 0.15, 0.175, 0.2, 0.25];  % Range of heights (m)
results.height_sweep.heights = height_range;
results.height_sweep.mass = zeros(size(height_range));
results.height_sweep.valid = zeros(size(height_range));

% Fixed parameters for height sweep
load_val = 9.81;   % Applied load (N·m)
theta_deg = 5;     % Pivot axis angle (degrees)
platform_diam = 0.14;  % Platform diameter (m)

for i = 1:length(height_range)
    fprintf('\nTesting Height = %.3f m\n', height_range(i));
    [mass, isValid] = Deflection_Optim(load_val, height_range(i), theta_deg, platform_diam, false);
    results.height_sweep.mass(i) = mass;
    results.height_sweep.valid(i) = isValid;
    fprintf('  Mass: %.2f g, Valid: %d\n', mass, isValid);
end

% Generate plots for the parameter sweeps
figure;
subplot(3,1,1);
plot(results.load_sweep.loads, results.load_sweep.mass, 'o-', 'LineWidth', 2);
xlabel('Load (N·m)');
ylabel('Optimized Mass (g)');
title('Effect of Load on Optimized Mass');
grid on;

subplot(3,1,2);
plot(results.angle_sweep.angles, results.angle_sweep.mass, 'o-', 'LineWidth', 2);
xlabel('Pivot Angle (degrees)');
ylabel('Optimized Mass (g)');
title('Effect of Pivot Angle on Optimized Mass');
grid on;

subplot(3,1,3);
plot(results.height_sweep.heights, results.height_sweep.mass, 'o-', 'LineWidth', 2);
xlabel('Mechanism Height (m)');
ylabel('Optimized Mass (g)');
title('Effect of Mechanism Height on Optimized Mass');
grid on;

% Plot validity regions
figure;
subplot(3,1,1);
bar(results.load_sweep.loads, results.load_sweep.valid, 0.6);
xlabel('Load (N·m)');
ylabel('Valid Design (1=Yes, 0=No)');
title('Valid Design Regions - Load Sweep');
ylim([0, 1.2]);
grid on;

subplot(3,1,2);
bar(results.angle_sweep.angles, results.angle_sweep.valid, 0.6);
xlabel('Pivot Angle (degrees)');
ylabel('Valid Design (1=Yes, 0=No)');
title('Valid Design Regions - Angle Sweep');
ylim([0, 1.2]);
grid on;

subplot(3,1,3);
bar(results.height_sweep.heights, results.height_sweep.valid, 0.6);
xlabel('Mechanism Height (m)');
ylabel('Valid Design (1=Yes, 0=No)');
title('Valid Design Regions - Height Sweep');
ylim([0, 1.2]);
grid on;

% Run a detailed analysis for a sample case with plots
fprintf('\n===== DETAILED ANALYSIS FOR SAMPLE CASE =====\n');
sample_load = 9.81;
sample_height = 0.15;
sample_angle = 5;
fprintf('Running detailed analysis for:\n');
fprintf('  Load: %.2f N·m\n', sample_load);
fprintf('  Height: %.3f m\n', sample_height);
fprintf('  Angle: %.2f degrees\n', sample_angle);
[mass, isValid] = Deflection_Optim(sample_load, sample_height, sample_angle, platform_diam, true);
fprintf('FINAL RESULTS:\n');
fprintf('  Optimized Mass: %.2f g\n', mass);
fprintf('  Valid Design: %d\n', isValid);

% Additional test: multi-parameter grid search for critical combinations
fprintf('\n===== GRID SEARCH FOR CRITICAL COMBINATIONS =====\n');
% Create a smaller grid to find boundaries
critical_loads = [5, 10, 15];
critical_angles = [5, 15, 30];

% Create a table to store results
fprintf('Load (N·m) | Angle (°) | Mass (g) | Valid\n');
fprintf('------------|-----------|----------|------\n');

for l = 1:length(critical_loads)
    for a = 1:length(critical_angles)
        [mass, isValid] = Deflection_Optim(critical_loads(l), H_val, critical_angles(a), platform_diam, false);
        fprintf('   %.1f    |    %.1f   |  %.1f  |  %d\n', critical_loads(l), critical_angles(a), mass, isValid);
    end
end

% Save the results structure
save('deflection_optim_results.mat', 'results');
fprintf('\nResults saved to deflection_optim_results.mat\n');