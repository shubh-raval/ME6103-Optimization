lb = [0.02, 0.01, 0]; % Lower bounds
ub = [.12, 1.5, 10000]; % Upper bounds
x0 = [.02, 1, 5]; % Initial guess

f1_max = 0.363480;
f2_max = 1.303012; 

weights = [0.1, .9; 0.3, 0.7; 0.5, 0.5; 0.7, 0.3; 0.9, 0.1];

% Store solutions
solutions = zeros(size(weights,1), 9); % Store w1, w2, D_p, L_e, N, Z, f1, f2,c

% Open file to save results
fid = fopen('weighted_solutions.txt', 'w');
fprintf(fid, 'w1\tw2\tD_p\tL_e\tN\tZ\tf1\tf2\tc\n');

% Loop through weight combinations and optimize
for i = 1:size(weights,1)
    w1 = weights(i,1);
    w2 = weights(i,2);

    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    
    weighted_obj = @(x) weighted_sum(x, f1_max, f2_max, w1, w2); 

    [x_opt, Z_opt] = fmincon(weighted_obj, x0, [], [], [], [], lb, ub, [], options); 

    [f1_opt, f2_opt, c] = compute_objectives(x_opt);

    solutions(i, :) = [w1, w2, x_opt(1), x_opt(2), x_opt(3), Z_opt, f1_opt, f2_opt, c];

    fprintf(fid, '%.1f\t%.1f\t%.4f\t%.4f\t%.4f\t%.6f\t%.6f\t%.6f\t%.6f\n', w1, w2, x_opt(1), x_opt(2), x_opt(3), Z_opt, f1_opt, f2_opt, c);
end


fclose(fid);


% Load pre-run data from solutions.txt
data = dlmread('solution.txt', '', 1, 0);
f1_values = data(:, 4);
f2_values = data(:, 5);

% Plot Pareto front
figure;
scatter(f1_values, f2_values, 'bo', 'DisplayName', 'Pre-run Solutions'); 
hold on;

% Plot optimized solutions from weighted sum
scatter(solutions(:, 7), solutions(:, 8), 100, 'r', 'filled', 'DisplayName', 'Weighted Sum Solutions');

% Labels and legend
xlabel('Objective 1 (f1)');
ylabel('Objective 2 (f2)');
title('Objective Trade-off: Pre-run vs. Weighted Sum Optimization');
legend;
grid on;
hold off;

function [f1, f2,c] = compute_objectives(x)
    D_p = x(1);
    L_e = x(2);
    N = x(3);

    % Dependent Variables
    S_p = D_p; 
    H_e = 1.1 * D_p;
    W_e = N * (2 * D_p) + D_p;
    Half_Circ = D_p*pi*.5;

    % Objective Function 1: Volume
    f1 = W_e * H_e * L_e;

    % Objective Function 2: Head Loss
    rho_water = 998; % kg/m^3
    m_dot = 0.2; % kg/s
    f_friction = 0.04;
    flow_v = m_dot / (rho_water * pi * (D_p / 2)^2);
    mhl_t1 = (flow_v^2) * (1 / (2 * 9.81));

    % Compute Straight Pipe Lengths
    ExposedStraightPipe_Len = (N + 1) * (L_e - 1.5 * D_p);
    if N > 1
        ExposedStraightPipe_Len = ExposedStraightPipe_Len + (N - 1) * (L_e - 3 * D_p);
    end
    U_bend_Len = N * (1.5 * D_p * pi * 0.5);

    % Total Head Loss Calculation
    major_head_loss = f_friction * (ExposedStraightPipe_Len / D_p) * mhl_t1;
    k = f_friction * 50; % Equivalent length ratio for return bends
    minor_head_loss = k * (flow_v^2) / (2 * 9.81);
    f2 = major_head_loss + N * minor_head_loss;
    SA = Half_Circ * (U_bend_Len + ExposedStraightPipe_Len);
    c = SA - 2.5104;
end


