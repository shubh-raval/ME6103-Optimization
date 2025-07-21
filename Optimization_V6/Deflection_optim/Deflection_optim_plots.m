% Surface Plot Script for Deflection Optimization
% This script calls Deflection_Optim and generates surface plots

% Define input parameters in the command window before running:
% load, H, theta_deg, platform_diam

% Check if the required inputs exist
if ~exist('load', 'var') || ~exist('H', 'var') || ~exist('theta_deg', 'var') || ~exist('platform_diam', 'var')
    error('Please define load, H, theta_deg, and platform_diam in the command window.');
end

% Run the optimization
[mass, isValid, W_opt, T_opt, a_opt, b_opt] = Deflection_Optim(load, H, theta_deg, platform_diam);

% Material properties
E = 69e9;          % Young's modulus (Pa)
density = 2700;    % Density in kg/m^3

% Calculate base radius and quarter ellipse length
base_rad = H * sind(theta_deg);
quarter_circumference = (pi/4) * (3*(a_opt + b_opt) - sqrt((3*a_opt + b_opt)*(a_opt + 3*b_opt)));

% Create mesh grid around optimal values
W_range = linspace(max(0.005, W_opt*0.5), min(0.05, W_opt*1.5), 30);
T_range = linspace(max(0.005, T_opt*0.5), min(0.05, T_opt*1.5), 30);
[W_grid, T_grid] = meshgrid(W_range, T_range);

% Calculate values across the grid
mass_grid = zeros(size(W_grid));
stiffness_grid = zeros(size(W_grid));
deflection_grid = zeros(size(W_grid));
valid_grid = false(size(W_grid));

for i = 1:size(W_grid, 1)
    for j = 1:size(W_grid, 2)
        W = W_grid(i, j);
        T = T_grid(i, j);
        
        % Calculate mass
        mass_grid(i, j) = density * W * T * quarter_circumference * 1000; % in grams
        
        % Moment of inertia
        I = (W * T^3) / 12;
        
        % Stiffness
        stiffness_grid(i, j) = E * I / (quarter_circumference^3);
        
        % Simplified deflection calculation
        R_eff = (H * (1 + cosd(theta_deg))) / 2;
        deflection_grid(i, j) = (load * R_eff * ((pi/2) - 1)) / (E * I); % in meters
        
        % Check if valid
        deflection_limit = 0.03 * quarter_circumference; % in meters
        valid_grid(i, j) = (deflection_grid(i, j) <= deflection_limit) && ...
                          (stiffness_grid(i, j) >= 2900) && (stiffness_grid(i, j) <= 20100);
    end
end

% Create figures
figure('Name', 'Optimization Results', 'Position', [100, 100, 1000, 800]);

% Mass surface
subplot(2, 2, 1);
surf(W_grid, T_grid, mass_grid);
hold on;
plot3(W_opt, T_opt, mass, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('Width (m)'); ylabel('Thickness (m)'); zlabel('Mass (g)');
title('Mass Surface');
colormap(jet); colorbar; grid on;
view(45, 30);

% Contour with valid region
subplot(2, 2, 2);
contourf(W_grid, T_grid, mass_grid, 20);
hold on;
contour(W_grid, T_grid, double(valid_grid), [0.5 0.5], 'r', 'LineWidth', 2);
plot(W_opt, T_opt, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('Width (m)'); ylabel('Thickness (m)');
title('Valid Design Region');
colormap(jet); colorbar; grid on;

% Stiffness surface
subplot(2, 2, 3);
surf(W_grid, T_grid, stiffness_grid);
hold on;
plot3(W_opt, T_opt, E * (W_opt * T_opt^3 / 12) / (quarter_circumference^3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('Width (m)'); ylabel('Thickness (m)'); zlabel('Stiffness (N/m)');
title('Stiffness Surface');
colormap(jet); colorbar; grid on;
view(45, 30);

% Deflection surface
subplot(2, 2, 4);
surf(W_grid, T_grid, deflection_grid);
hold on;
R_eff_opt = (H * (1 + cosd(theta_deg))) / 2;
I_opt = (W_opt * T_opt^3) / 12;
deflection_opt = (load * R_eff_opt * ((pi/2) - 1)) / (E * I_opt);
plot3(W_opt, T_opt, deflection_opt, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('Width (m)'); ylabel('Thickness (m)'); zlabel('Deflection (m)');
title('Deflection Surface');
colormap(jet); colorbar; grid on;
view(45, 30);

% Print results
fprintf('\n======= SURFACE PLOT RESULTS =======\n');
fprintf('Optimized values:\n');
fprintf('Mass: %.2f grams\n', mass);
fprintf('W: %.5f m\n', W_opt);
fprintf('T: %.5f m\n', T_opt);
fprintf('a: %.5f m\n', a_opt);
fprintf('b: %.5f m\n', b_opt);
fprintf('===================================\n\n');