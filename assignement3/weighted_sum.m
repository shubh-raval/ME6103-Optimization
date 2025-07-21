function Z = weighted_sum(x, f1_max, f2_max,w1,w2)
    % Extract decision variables
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
    Total_Len =   U_bend_Len + ExposedStraightPipe_Len; % m
    Total_surface_area = Total_Len * Half_Circ;
    required_SA = 2.5104;
    % Total Head Loss Calculation
    major_head_loss = f_friction * (ExposedStraightPipe_Len / D_p) * mhl_t1;
    k = f_friction * 50; % Equivalent length ratio for return bends
    minor_head_loss = k * (flow_v^2) / (2 * 9.81);
    f2 = major_head_loss + N * minor_head_loss; % Objective 2: Head Loss
    % If constraint is violated, add penalty
    if Total_surface_area < 2.5104
        penalty = 1e6 * (required_SA - Total_surface_area); 
        f1 = f1 + penalty;
        f2 = f2 + penalty;
    end
    % If constraint is overly exceeded, add penalty
    if Total_surface_area > 2.5104 * 1.1
        penalty = 1e6 * (Total_surface_area - required_SA);
        f1 = f1 + penalty;
        f2 = f2 + penalty;
    end


    Z = (w1*f1 / f1_max) + (w2*f2 / f2_max);

end
