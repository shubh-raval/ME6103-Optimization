% test_inverse_kinematics.m
clear; clc;

% Define test parameters:
params.H = 0.15;             % Mechanism height (m)
base_rad = 0.12;             % Base platform radius (m)
alpha1 = atand(base_rad / params.H);  % Compute alpha1 (should be ~38.66°)
alpha2 = 50;               % Set alpha2 to 20° for testing
beta2 = 90;                % Set beta2 to 90° (as used in your kinematics functions)
eta11 = 0; eta12 = 120; eta13 = 240;  % Base joint offset angles

% Define a test orientation:
test_phi = 10;
test_theta = 5;
test_psi = 0;

% Call calculate_vi with the test orientation for leg 1:
v = calculate_vi(test_phi, test_theta, test_psi, beta2, eta11, eta12, eta13, 1);

% Compute the coefficients for the first leg using define_coefficients:
[A, B, C] = define_coefficients(alpha1, eta11, eta12, eta13, alpha2, v, 1);

% Calculate the discriminant for the quadratic equation:
discriminant = ((2 * B)^2) - 4 * A * C;
fprintf('Discriminant = %.4f\n', discriminant);

if discriminant < 0
    error('Discriminant is negative, inverse kinematics cannot be solved for this configuration.');
end

% Solve for t1 using the tangent-half-angle method:
t1 = (-(2 * B) + sqrt(discriminant)) / (2 * A);

% Compute the joint angle for leg 1 (in degrees):
test_theta_i = rad2deg(2 * atan2(t1, 1));

disp(['Test inverse kinematics solution for leg 1: ', num2str(test_theta_i), ' degrees']);