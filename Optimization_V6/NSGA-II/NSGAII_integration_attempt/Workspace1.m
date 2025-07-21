function [W, N, nn, CN] = Workspace1(alpha1, eta11, eta12, eta13, alpha2, beta2, gamma2, phi, theta, psi, zeta_min)
% Workspace1 computes the valid workspace configurations and conditioning
% indices for the SPM.
%
% Outputs:
%   W  - Matrix of valid Euler angle configurations.
%   N  - Matrix of corresponding end-effector z-axis directions.
%   nn - Matrix of singular configurations (if any).
%   CN - Vector of conditioning indices at each valid configuration.
%
% This function uses inverse kinematics and the Jacobian to determine the
% conditioning index for each configuration.

W = [];
N = [];
nn = [];
CN = [];

% Debug info
fprintf('Starting Workspace1 calculation:\n');
fprintf('  alpha1=%.2f, alpha2=%.2f, beta2=%.2f, gamma2=%.2f\n', alpha1, alpha2, beta2, gamma2);
fprintf('  Grid size: phi=%d points, theta=%d points, psi=%d points\n', length(phi), length(theta), length(psi));
fprintf('  zeta_min=%.4f\n', zeta_min);

% Make sure zeta_min is not too restrictive
if zeta_min > 0.1
    fprintf('WARNING: zeta_min (%.4f) may be too restrictive. Consider values ≤ 0.05\n', zeta_min);
end

% Count processed configurations for debugging
total_configs = 0;
valid_configs = 0;
singular_configs = 0;
failed_ik = 0;

for i = 1:length(phi)
    for j = 1:length(theta)
        for k = 1:length(psi)
            total_configs = total_configs + 1;
            
            % Try to compute inverse kinematics solution
            try
                T = inverse_kinematics(alpha1, eta11, eta12, eta13, alpha2, beta2, phi(i), theta(j), psi(k));
                
                % Skip if any joint angle is NaN or Inf
                if any(isnan(T)) || any(isinf(T))
                    failed_ik = failed_ik + 1;
                    continue;
                end
                
                % We have a valid IK solution
                valid_configs = valid_configs + 1;
                
                % Store Euler angles and compute platform orientation
                w = [phi(i); theta(j); psi(k)];
                Q = define_Q(phi(i), theta(j), psi(k));
                n = Q(:,3);
                W = [W, w];
                N = [N, n];
                
                % Compute forward kinematics and Jacobian
                [v1, v2, v3, w1, w2, w3] = forward_kinematics(alpha1, eta11, eta12, eta13, alpha2, gamma2, T(1), T(2), T(3));
                
                % Check for potential numerical issues in vectors
                if any(isnan([v1; v2; v3; w1; w2; w3])) || any(isinf([v1; v2; v3; w1; w2; w3]))
                    failed_ik = failed_ik + 1;
                    continue;
                end
                
                % Calculate Jacobian and conditioning index
                try
                    J = calculate_Jacobian(w1, w2, w3, v1, v2, v3);
                    zeta = conditioning_index(J);
                    
                    % Handle NaN or negative conditioning index
                    if isnan(zeta) || ~isreal(zeta) || zeta < 0
                        singular_configs = singular_configs + 1;
                        nn = [nn, n];
                        continue;
                    end
                    
                    % Add conditioning index to list
                    CN = [CN, zeta];
                    
                    % Check if this is a near-singular configuration
                    if zeta < zeta_min
                        singular_configs = singular_configs + 1;
                        nn = [nn, n];
                    end
                catch e
                    % Handle errors in Jacobian calculation
                    failed_ik = failed_ik + 1;
                    fprintf('Error calculating Jacobian for phi=%f, theta=%f, psi=%f: %s\n', ...
                           phi(i), theta(j), psi(k), e.message);
                end
            catch e
                % Handle errors in inverse kinematics
                failed_ik = failed_ik + 1;
                if mod(failed_ik, 100) == 1  % Limit output to avoid flooding console
                    fprintf('Error in inverse kinematics for phi=%f, theta=%f, psi=%f: %s\n', ...
                           phi(i), theta(j), psi(k), e.message);
                end
            end
        end
    end
end

% Print comprehensive debug information
fprintf('\nWorkspace Analysis Summary:\n');
fprintf('  Total configurations processed: %d\n', total_configs);
fprintf('  Valid IK solutions found: %d (%.1f%%)\n', valid_configs, 100*valid_configs/total_configs);
fprintf('  Failed IK or FK calculations: %d (%.1f%%)\n', failed_ik, 100*failed_ik/total_configs);
fprintf('  Near-singular configurations: %d\n', singular_configs);
fprintf('  Valid conditioning indices found: %d\n', length(CN));

% Additional stats on conditioning indices if available
if ~isempty(CN)
    fprintf('  Conditioning index statistics:\n');
    fprintf('    Minimum: %.6f\n', min(CN));
    fprintf('    Maximum: %.6f\n', max(CN));
    fprintf('    Mean: %.6f\n', mean(CN));
    fprintf('    Median: %.6f\n', median(CN));
else
    fprintf('WARNING: No valid conditioning indices found! The mechanism has no valid workspace for these parameters.\n');
    fprintf('Try reducing zeta_min or adjusting geometric parameters.\n');
end

end