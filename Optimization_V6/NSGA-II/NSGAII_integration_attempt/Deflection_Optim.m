function [mass, isValid, W_opt, T_opt] = Deflection_Optim(load, h, theta_deg, platform_diam, debug, W_init, T_init)
% DEFLECTION_OPTIM Optimizes the SPM linkage design by minimizing mass subject to 
% deflection and stiffness constraints under a specified torque load.
%
%   [mass, isValid, W_opt, T_opt] = Deflection_Optim(load, h, theta_deg, platform_diam, debug, W_init, T_init)
%
%   Inputs:
%       load         - Applied load (torque in N·m) from Simscape.
%       h            - Mechanism height (m); used as the major half-axis (a) of the ellipse.
%       theta_deg    - Pivot axis angle (deg); used to compute the minor half-axis (b = h*cosd(theta_deg)).
%       platform_diam- Platform diameter (m).
%       debug        - Boolean flag; if true, prints debug info and generates plots.
%       W_init       - Initial guess for link width (m). If not provided, a random value is used.
%       T_init       - Initial guess for link thickness (m). If not provided, a random value is used.
%
%   Outputs:
%       mass         - Optimized mass (in grams) computed as:
%                      mass = density * (W * T * L) * 1000,
%                      where L is the quarter ellipse arc length.
%       isValid      - Boolean flag: true if the computed deflection is within limits and 
%                      stiffness is within the specified bounds; false otherwise.
%       W_opt        - Optimized link width (m)
%       T_opt        - Optimized link thickness (m)
%
%   Notes:
%       - Ellipse half-axes:
%             a = h  (major half-axis)
%             b = h * cosd(theta_deg)  (minor half-axis)
%
%       - Quarter ellipse arc length (L) is computed using Ramanujan's approximation.
%       - Effective radius for curved beam: R_eff = (h*(1+cosd(theta_deg)))/2.
%       - Deflection: deflection = (M * R_eff / (E * I))*(pi/2 - 1), where I = W*T^3/12.
%       - Constraints:
%             deflection <= 0.03 * L
%             stiffness K = E*I/L^3 must be between K_min and K_max.
%
%       - Material properties (example for aluminum):
%             E = 69e9 (Pa)
%             density = 2700 (kg/m^3)
%
%       - An adaptive penalty factor is applied if constraints are violated.
%
%       - debug: Controls detailed output and plots.

% Material properties
E = 69e9;          % Young's modulus in Pascals
density = 2700;    % Density in kg/m^3

% --- Compute Ellipse Half-Axes from h and theta_deg ---
a = h;                   % Major half-axis (m)
b = h * cosd(theta_deg); % Minor half-axis (m)
% ----------------------------------------------------------

% Compute the quarter ellipse arc length using Ramanujan's approximation
L = approximate_quarter_ellipse_length(a, b);

% Maximum allowable deflection: 3% of the effective linkage length L
deflection_limit = 0.03 * L;

% --- Define Torque Load ---
M = load;  % (N·m)

% --- Compute Effective Radius for the Curved Beam ---
R_eff = (h * (1 + cosd(theta_deg))) / 2;

% --- Set Bounds for Design Variables (W and T) ---
min_val = 0.005;  % Lower bound in meters
max_val = 0.5;   % Upper bound in meters

% --- Initialize Design Guesses for W and T ---
if nargin < 6 || isempty(W_init)
    W = min_val + (max_val - min_val) * rand();
else
    W = W_init;
end

if nargin < 7 || isempty(T_init)
    T = min_val + (max_val - min_val) * rand();
else
    T = T_init;
end

% --------------------------------------------------------------

% Set optimization parameters
max_iterations = 20;  % Reduced maximum iterations for gradient descent
init_learning_rate = 1e-4;
learning_rate = init_learning_rate;
min_learning_rate = 1e-7;
learning_rate_decay = 0.95;

% Convergence criteria
conv_tolerance = 1e-6;  
consecutive_stable_iterations = 0;  
required_stable_iterations = 5;  

prev_mass_vol = Inf;

% Define a tolerance for constraint comparisons (adjusted to 1e-6)
tol = 1e-6;

% Define a flexible penalty factor (scaling factor) that controls how steeply the penalty increases.
% A higher value makes the penalty steeper. Here we choose 1e3 as a moderate value.
penalty_factor = 1e3;

% (Optional) History storage for debug plotting
if debug
    WHistory = zeros(max_iterations, 1);
    THistory = zeros(max_iterations, 1);
    massHistory = zeros(max_iterations, 1);
    deflHistory = zeros(max_iterations, 1);
    stiffHistory = zeros(max_iterations, 1);
    learningRateHistory = zeros(max_iterations, 1);
end

% Gradient descent loop: adjust W and T to minimize mass while meeting constraints.
for i = 1:max_iterations
    % Compute current mass surrogate (volume in m^3)
    current_mass_vol = W * T * L;
    
    % Compute moment of inertia (I = W*T^3/12)
    I = W * T^3 / 12;
    
    % Compute deflection using the curved beam formula:
    deflection_est = (M * R_eff * ((pi/2) - 1)) / (E * I);
    
    % Compute stiffness K = E*I / L^3 (N/m)
    K = E * I / (L^3);
    K_min = 3000;   % Minimum required stiffness (N/m)
    K_max = 20000;  % Upper stiffness limit (N/m)
    
    % Save iteration history if debugging is enabled
    if debug
        WHistory(i) = W;
        THistory(i) = T;
        massHistory(i) = current_mass_vol;
        deflHistory(i) = deflection_est;
        stiffHistory(i) = K;
        learningRateHistory(i) = learning_rate;
    end
    
    % Check convergence based on relative change in mass surrogate
    rel_change = abs(current_mass_vol - prev_mass_vol) / max(prev_mass_vol, eps);
    if rel_change < conv_tolerance
        consecutive_stable_iterations = consecutive_stable_iterations + 1;
        if consecutive_stable_iterations >= required_stable_iterations
            if debug
                fprintf('Converged after %d iterations (relative change < %e for %d consecutive iterations)\n', i, conv_tolerance, required_stable_iterations);
            end
            break;
        end
    else
        consecutive_stable_iterations = 0;
    end
    prev_mass_vol = current_mass_vol;
    
    % Compute constraint violations
    deflection_violation = max(0, deflection_est - deflection_limit);
    stiffness_violation = max(0, K_min - K);  % violation if too low
    upper_stiffness_violation = max(0, K - K_max);  % violation if too high
    
    % Calculate adaptive penalty factors using the penalty_factor:
    penalty_deflection = 1 + penalty_factor * (deflection_violation / deflection_limit);
    penalty_stiffness_low = 1 + penalty_factor * (stiffness_violation / K_min);
    penalty_stiffness_high = 1 + penalty_factor * (upper_stiffness_violation / K_max);
    
    % Update design variables:
    if deflection_violation > 0 || stiffness_violation > 0 || upper_stiffness_violation > 0
        if upper_stiffness_violation > 0
            % If stiffness is too high, reduce dimensions
            W_new = W - learning_rate * W * 0.5 * penalty_stiffness_high;
            T_new = T - learning_rate * T * 1.0 * penalty_stiffness_high;
        else
            % If deflection is too high or stiffness is too low, increase dimensions
            W_new = W + learning_rate * W * 0.5 * (penalty_deflection + penalty_stiffness_low);
            T_new = T + learning_rate * T * 1.0 * (penalty_deflection + penalty_stiffness_low);
        end
    else
        % Constraints satisfied: try to reduce mass
        W_new = W - learning_rate * W * 1.2;
        T_new = T - learning_rate * T * 0.8;
    end
    
    % Apply bounds to ensure values are within [min_val, max_val]
    W_new = max(min_val, min(max_val, W_new));
    T_new = max(min_val, min(max_val, T_new));
    
    % Adapt learning rate based on progress
    if i > 1 && (deflection_violation > 0 || stiffness_violation > 0 || upper_stiffness_violation > 0)
        learning_rate = max(min_learning_rate, learning_rate * learning_rate_decay);
    elseif mod(i, 10) == 0
        learning_rate = max(min_learning_rate, learning_rate * learning_rate_decay);
    end
    
    W = W_new;
    T = T_new;
end

% Compute final mass surrogate and final constraints
final_mass_vol = W * T * L;
I = W * T^3 / 12;
final_deflection = (M * R_eff * ((pi/2) - 1)) / (E * I);
K = E * I / (L^3);

% --- Adaptive Penalty Application for Deflection with tolerance ---
if final_deflection > deflection_limit + tol
    penalty_deflection_final = 1 + penalty_factor * ((final_deflection - deflection_limit) / deflection_limit);
    final_mass_vol = final_mass_vol * penalty_deflection_final;
    if debug
        fprintf('Deflection constraint not met: %.6f m (Limit: %.6f m). Penalty factor: %.2f\n', final_deflection, deflection_limit, penalty_deflection_final);
    end
    isValid = false;
else
    isValid = true;
end

% --- Adaptive Penalty Application for Stiffness (Lower Bound) ---
if K < K_min - tol
    penalty_stiffness_low_final = 1 + penalty_factor * ((K_min - K) / K_min);
    final_mass_vol = final_mass_vol * penalty_stiffness_low_final;
    if debug
        fprintf('Stiffness constraint not met (low): K = %.2f N/m (Required: >= %.2f N/m). Penalty factor: %.2f\n', K, K_min, penalty_stiffness_low_final);
    end
    isValid = false;
end

% --- Adaptive Penalty Application for Stiffness (Upper Bound) ---
if K > K_max + tol
    penalty_stiffness_high_final = 1 + penalty_factor * ((K - K_max) / K_max);
    final_mass_vol = final_mass_vol * penalty_stiffness_high_final;
    if debug
        fprintf('Upper stiffness constraint not met: K = %.2f N/m (Maximum allowed: %.2f N/m). Penalty factor: %.2f\n', K, K_max, penalty_stiffness_high_final);
    end
    isValid = false;
end

% Convert volume to mass (in grams)
mass = density * final_mass_vol * 1000;

% Return optimized W and T
W_opt = W;
T_opt = T;

% --- Debug Prints and Plots ---
if debug
    fprintf('Final design parameters:\n');
    fprintf('   W (width)     = %.5f m\n', W);
    fprintf('   T (thickness) = %.5f m\n', T);
    fprintf('Optimized mass surrogate (W*T*L, in m^3): %.6e m^3\n', final_mass_vol);
    fprintf('Final deflection: %.6f m (Limit: %.6f m)\n', final_deflection, deflection_limit);
    fprintf('Stiffness: %.2f N/m (Required: [%.2f, %.2f] N/m)\n', K, K_min, K_max);
    fprintf('Optimized Mass (in grams): %.6f\n', mass);
    
    actual_iterations = min(max_iterations, length(WHistory));
    
    figure;
    subplot(5,1,1);
    plot(WHistory(1:actual_iterations), 'LineWidth', 2);
    xlabel('Iteration'); ylabel('W (m)');
    title('Width (W) Evolution'); grid on;
    
    subplot(5,1,2);
    plot(THistory(1:actual_iterations), 'LineWidth', 2);
    xlabel('Iteration'); ylabel('T (m)');
    title('Thickness (T) Evolution'); grid on;
    
    subplot(5,1,3);
    plot(massHistory(1:actual_iterations), 'LineWidth', 2);
    xlabel('Iteration'); ylabel('Mass Surrogate (m^3)');
    title('Mass Surrogate Evolution'); grid on;
    
    subplot(5,1,4);
    plot(deflHistory(1:actual_iterations), 'LineWidth', 2);
    hold on; plot([1, actual_iterations], [deflection_limit, deflection_limit], 'r--', 'LineWidth', 1.5);
    legend('Deflection','Limit');
    xlabel('Iteration'); ylabel('Deflection (m)');
    title('Deflection Evolution'); grid on;
    
    subplot(5,1,5);
    semilogy(learningRateHistory(1:actual_iterations), 'LineWidth', 2);
    xlabel('Iteration'); ylabel('Learning Rate');
    title('Learning Rate Adaptation'); grid on;
    
    figure(1);
    W_min_plot = min(WHistory(1:actual_iterations)) * 0.95;
    W_max_plot = max(WHistory(1:actual_iterations)) * 1.05;
    T_min_plot = min(THistory(1:actual_iterations)) * 0.95;
    T_max_plot = max(THistory(1:actual_iterations)) * 1.05;
    
    [W_grid, T_grid] = meshgrid(linspace(W_min_plot, W_max_plot, 50), linspace(T_min_plot, T_max_plot, 50));
    mass_grid = W_grid .* T_grid * L;
    
    figure(2);
    surf(W_grid, T_grid, mass_grid, 'FaceAlpha', 0.7);
    hold on;
    plot3(WHistory(1:actual_iterations), THistory(1:actual_iterations), massHistory(1:actual_iterations), 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    xlabel('W (m)'); ylabel('T (m)'); zlabel('Mass Surrogate (m^3)');
    title('Gradient Descent Trajectory: Evolution of W and T'); grid on;
    view(135,30);
    
    figure(3);
    subplot(2,1,1);
    plot(stiffHistory(1:actual_iterations), 'LineWidth', 2);
    hold on; plot([1, actual_iterations], [K_min, K_min], 'r--', 'LineWidth', 1.5);
    legend('Stiffness','Minimum Required');
    xlabel('Iteration'); ylabel('Stiffness (N/m)');
    title('Stiffness Evolution'); grid on;
    
    subplot(2,1,2);
    yyaxis left;
    plot(massHistory(1:actual_iterations), 'b-', 'LineWidth', 2);
    ylabel('Mass Surrogate (m^3)');
    yyaxis right;
    plot(deflHistory(1:actual_iterations)./deflection_limit, 'r-', 'LineWidth', 2);
    hold on; plot([1, actual_iterations], [1, 1], 'r--', 'LineWidth', 1.5);
    ylabel('Deflection Ratio (actual/limit)');
    xlabel('Iteration');
    title('Convergence Progress'); grid on;
end

end

function arc_len = approximate_quarter_ellipse_length(a, b)
% APPROXIMATE_QUARTER_ELLIPSE_LENGTH Computes the quarter perimeter of an ellipse
% using Ramanujan's approximation.
%
%   arc_len = approximate_quarter_ellipse_length(a, b)
%
%   Inputs:
%       a - Major half-axis (m)
%       b - Minor half-axis (m)
%
%   Ramanujan's approximation for the full ellipse perimeter is:
%       P ≈ π * (3*(a+b) - sqrt((3*a+b)*(a+3*b)))
%   The quarter ellipse arc length is:
%       arc_len = P / 4;
%
%   Units: meters.
    perimeter = pi * (3*(a + b) - sqrt((3*a + b) * (a + 3*b)));
    arc_len = perimeter / 4;
end
