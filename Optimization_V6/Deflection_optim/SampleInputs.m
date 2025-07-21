% Sample test parameters from Simscape
load_val = 9.81;      % Applied load (torque in N·m); e.g., for a 5 kg payload with 2× factor
H_val = 0.15;         % Mechanism height (m) (major half-axis)
theta_deg = 5;        % Pivot axis angle in degrees
platform_diam = 0.14; % Platform diameter in meters (e.g., 140 mm)
debug = true;         % Enable debug mode to display plots and iteration details

% Call the optimization function
[mass, isValid] = Deflection_Optim(load_val, H_val, theta_deg, platform_diam, debug);

% Display the results
fprintf('Optimized Mass: %.6f grams\n', mass);
fprintf('Design is Valid (deflection & stiffness within permissible range): %d\n', isValid);
