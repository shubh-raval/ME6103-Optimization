% Sample input values:
load_val = 9.81; % Applied load (N·m) 
h_val = 0.15; % Mechanism height (m)
theta_deg_val = 5; % Pivot axis angle (deg) 
platform_diam_val = 0.14; % Platform diameter (m)
debug_flag = false; % Debug flag (set to true for more info/plots)
W_init = 0.01; % Initial guess for link width (m)
T_init = 0.01; % Initial guess for link thickness (m)

% Run deflection optimization:
[mass, isValid, W_opt, T_opt] = Deflection_Optim(load_val, h_val, theta_deg_val, platform_diam_val, debug_flag, W_init, T_init);
fprintf('Optimized mass: %.2f grams\n', mass);
fprintf('Optimized W: %.5f m, Optimized T: %.5f m\n', W_opt, T_opt);

% Define additional parameters for kinematics:
base_rad = 1.2; % Base platform radius (m)
top_rad = platform_diam_val / 2; % Mobile platform radius (m) 
h = h_val; % Use h_val for mechanism height
theta = theta_deg_val; % Pivot axis angle (deg)
alpha1 = atand(base_rad / h); % Compute alpha1 (deg)
alpha2 = []; % Let SPM_kinematics compute alpha2 internally

% Call SPM_kinematics and get data for plotting:
[GCI, x_data, y_data, z_data, CN] = SPM_kinematics(base_rad, top_rad, h, theta, alpha1, alpha2, isValid);
fprintf('Calculated GCI: %.4f\n', GCI);

% Plot workspace
figure;
plot3(x_data, y_data, z_data, '.');
hold on;
xlabel('x'); ylabel('y'); zlabel('z');
zlim([0,1]); grid on; grid minor; set(gcf, 'color', 'w');
legend('non-singular');
str_title = sprintf('Range of Motion: h = %.2f m, \\theta = %d\\circ', h, theta);
title(str_title);

% Plot conditioning number
figure;
colorpost = CN';
scatter3(x_data, y_data, z_data, 5 * ones(size(x_data)), colorpost, 'filled');
c = colorbar;
ylabel(c, 'Conditioning #', 'FontSize', 8, 'Rotation', 270);
xlim([-1 1]); ylim([-1 1]); zlim([0 1]); grid on; grid minor; set(gcf, 'color', 'w');
xlabel('x'); ylabel('y'); zlabel('z'); clim([0, 1]);
str_title2 = sprintf('Range of Motion w/ GCI: h = %.2f m, \\theta = %d\\circ', h, theta);
title(str_title2);