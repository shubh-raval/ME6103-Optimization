function GCI = SPM_kinematics_silent(base_rad, top_rad, h, theta_input, alpha1, alpha2, isValid)
    % Display input parameters for debugging
    fprintf('SPM_kinematics_silent inputs: base_rad=%.4f, top_rad=%.4f, h=%.4f, theta=%.4f\n', ...
            base_rad, top_rad, h, theta_input);
    
    % Fixed angular parameters for platform joints
    eta11 = 0; eta12 = 120; eta13 = 240;
    
    % IMPORTANT: Use the same beta2 and gamma2 formulation as in SPM_kinematics.m
    % to maintain consistency
    beta2 = 90;  % Fixed at 90 degrees
    gamma2 = 120;  % Fixed at 120 degrees
    
    % Calculate alpha1 if not provided (using h and base_rad)
    if isempty(alpha1) && ~isempty(h)
        alpha1 = atand(base_rad/h);
        fprintf('Calculated alpha1 = %.4f degrees\n', alpha1);
    end
    
    % Calculate alpha2 if not provided
    if isempty(alpha2) && ~isempty(theta_input)
        a_mag = sqrt(h^2 + base_rad^2);
        
        Ry_clockwise = [cosd(alpha1), 0, sind(alpha1);
                         0,           1,           0;
                        -sind(alpha1),0,  cosd(alpha1)];
        
        Rz_counter = [cosd(theta_input), -sind(theta_input), 0;
                       sind(theta_input),  cosd(theta_input), 0;
                       0,                  0,                 1];
        
        a_vec = Ry_clockwise * [0; 0; -1];
        a_vec = a_vec * a_mag;
        
        x_vec = [1; 0; 0];
        
        b_vec = Rz_counter * x_vec * top_rad;
        
        alpha2 = rad2deg(acos(dot(a_vec, b_vec) / (norm(a_vec) * norm(b_vec))));
        if alpha2 < 0
            alpha2 = alpha2 + 360;
        end
        fprintf('Calculated alpha2 = %.4f degrees\n', alpha2);
    end
    
    % Use a coarser grid for initial testing, but ensure adequate sampling
    phi = 0:40:360;     
    theta = 0:40:360;   
    psi = 0:40:360;    
    
    % Use very permissive zeta_min to find any valid configurations
    zeta_min = 0.01;    % Reduced from 0.05
    
    fprintf('Calculating workspace with zeta_min = %.4f (very permissive)\n', zeta_min);
    [W,N,nn,CN] = Workspace1(alpha1, eta11, eta12, eta13, alpha2, beta2, gamma2, phi, theta, psi, zeta_min);
    
    fprintf('Workspace1 returned %d configurations with conditioning indices\n', length(CN));
    
    % If no valid configurations found with permissive zeta_min, retry with even more permissive settings
    if isempty(CN)
        fprintf('No valid configurations found with zeta_min=%.4f, trying with zeta_min=0.001\n', zeta_min);
        zeta_min = 0.001;  % Extremely permissive
        
        % Try with a finer grid
        phi = 0:20:360;
        theta = 0:20:360;
        psi = 0:20:360;
        
        [W,N,nn,CN] = Workspace1(alpha1, eta11, eta12, eta13, alpha2, beta2, gamma2, phi, theta, psi, zeta_min);
        fprintf('Second attempt returned %d configurations\n', length(CN));
    end
    
    % If still no valid configurations, assign a small non-zero value
    if isempty(CN)
        GCI = 0.001;  % Small but non-zero
        fprintf('WARNING: No valid configurations found, setting GCI to %.4f\n', GCI);
    else
        % If at least one valid configuration found, compute mean GCI
        GCI = mean(CN);
        fprintf('Calculated mean GCI = %.6f from %d valid configurations\n', GCI, length(CN));
    end
    
    % Apply a less severe penalty for invalid deflection/stiffness designs
    if ~isValid
        old_GCI = GCI;
        GCI = GCI * 0.5;  % 50% penalty factor
        fprintf('Applied penalty factor 0.5 for invalid link design: %.6f -> %.6f\n', old_GCI, GCI);
    end
    
    % Ensure GCI is never exactly zero or negative
    if GCI <= 1e-6
        GCI = 1e-6;
        fprintf('GCI was too small or negative, set to minimum value: %.6f\n', GCI);
    end
    
    fprintf('Final GCI value: %.6f\n', GCI);
end