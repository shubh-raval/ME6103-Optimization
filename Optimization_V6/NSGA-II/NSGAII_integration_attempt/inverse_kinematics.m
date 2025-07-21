function [theta_i] = inverse_kinematics(alpha1, eta11, eta12, eta13, alpha2, beta2, phi, theta, psi)
% inverse_kinematics computes the input joint angles for a coaxial SPM.
%
% Inputs:
%   alpha1   - Base joint tilt angle (scalar)
%   eta11, eta12, eta13 - Base joint offset angles (scalars)
%   alpha2   - Platform joint tilt angle (scalar)
%   beta2    - Platform geometric parameter (scalar)
%   phi, theta, psi - Euler angles for platform orientation (scalars)
%
% Output:
%   theta_i  - 3x1 vector of computed joint angles (in degrees) for each leg.
%
% Note: This function computes one solution per leg.

% Initialize output as a 3x1 vector.
theta_i = zeros(3, 1);

for i = 1:3
    % Calculate the vector for the i-th leg given the platform orientation.
    v = calculate_vi(phi, theta, psi, beta2, eta11, eta12, eta13, i);
    
    % Compute the coefficients for the quadratic equation used in inverse kinematics.
    [A, B, C] = define_coefficients(alpha1, eta11, eta12, eta13, alpha2, v, i);
    
    % Handle numerical issues - ensure coefficients aren't too small
    if abs(A) < 1e-10
        % A is too close to zero - can cause numerical issues
        theta_i(i) = nan;
        continue;
    end
    
    % Compute the discriminant with robustness against numerical errors
    discriminant = (2*B)^2 - 4*A*C;
    
    % Numerical tolerance for discriminant - handle almost-zero cases
    if abs(discriminant) < 1e-10
        discriminant = 0;
    end
    
    % Check for negative discriminant (no real solutions)
    if discriminant < 0
        % No real solutions exist
        theta_i(i) = nan;
        continue;
    end
    
    % Compute tangent half-angle solution
    t_val = -2*B + sqrt(discriminant);
    t1 = t_val / (2*A);
    
    % Convert to joint angle, handling potential numerical issues
    if isnan(t1) || isinf(t1) || ~isreal(t1)
        theta_i(i) = nan;
    else
        theta_i(i) = rad2deg(2 * atan2(t1, 1));
    end
end

% Ensure joint angles are in a reasonable range (-180 to 180)
for i = 1:3
    if ~isnan(theta_i(i))
        % Normalize angle to -180 to 180 range
        theta_i(i) = mod(theta_i(i) + 180, 360) - 180;
    end
end
end