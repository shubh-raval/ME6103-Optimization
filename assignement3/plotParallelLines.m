function plotParallelLines(N, L_e, D_p)
    % Define the start and end x-coordinates for the lines
    x_start = 1.5*D_p;
    x_end = L_e - 1.5 * D_p;
    x_vals = [x_start, x_end];
    
    figure; hold on; grid on;
    title('Parallel Line Sets with Connecting Bends');
    xlabel('X-axis');
    ylabel('Y-axis');
    axis equal;
    
    % Define colors for different sets
    colors = ['b', 'r', 'g', 'm', 'c', 'k'];
    
    % Plot 2*(N+1) lines, as each set has 2 parallel lines
    for i = 0:N
        % Compute the y-coordinates for the set
        y1 = i * 2 * D_p;  % Offset each set correctly for 2 lines per set
        y2 = y1 + D_p;     % The second line in the set
        
        % Choose color from the list, cycling if needed
        color_index = mod(i, length(colors)) + 1;
        
        % Plot both lines in the set
        plot(x_vals, [y1, y1], 'Color', colors(color_index), 'LineWidth', 2);
        plot(x_vals, [y2, y2], 'Color', colors(color_index), 'LineWidth', 2);
    end
    
    % Keep hold on for additional plotting
    hold on;
    
    % Plot inner half-circle bends (centered at 1.5 * D_p) for even indices
    for i = 0:2:N-1  % Even indices only
        r_inner = 0.5 * D_p;
        y_inner_center = 1.5 * D_p + (i * 2 * D_p);  % Adjusted to start at 1.5*D_p
        theta = linspace(-pi/2, pi/2, 50);  % Half-circle (from -pi/2 to pi/2)
        
        % Calculate the x-coordinate of the inner bend based on parallel lines
        x_bend = x_end + 1.5*D_p;  % Set the x-coordinate after the last parallel line
        
        % Tangent to x = 1 for even i
        x_inner = x_bend + r_inner * cos(theta);  % Use standard half-circle
        
        % Adjust x_inner based on the parity of i
        x_inner = x_inner - 1.5 * D_p;  % Add 1.5 * D_p
        
        y_inner = y_inner_center + r_inner * sin(theta);
        plot(x_inner, y_inner, 'k-', 'LineWidth', 2);
    end
    
    % Plot inner half-circle bends (centered at 1.5 * D_p) for odd indices
    for i = 1:2:N-1  % Odd indices only
        r_inner = 0.5 * D_p;
        y_inner_center = 1.5 * D_p + (i * 2 * D_p);  % Adjusted to start at 1.5*D_p
        theta = linspace(-pi/2, pi/2, 50);  % Half-circle (from -pi/2 to pi/2)
        
        % Calculate the x-coordinate of the inner bend based on parallel lines
        x_bend = x_end + 1.5*D_p;  % Set the x-coordinate after the last parallel line
        
        % Tangent to x = 0 for odd i
        x_inner = x_bend - r_inner * cos(theta);  % Mirror half-circle (reverse direction)
        
        % Adjust x_inner based on the parity of i
        x_inner = x_inner + 1.5 * D_p - L_e ;  % Subtract 1.5 * D_p

        y_inner = y_inner_center + r_inner * sin(theta);
        plot(x_inner, y_inner, 'k-', 'LineWidth', 2);
    end
    
    % Plot outer half-circle bends (centered at 1.5 * D_p) for even indices
    for i = 0:2:N-1  % Even indices only
        r_outer = 1.5 * D_p;
        y_outer_center = 1.5 * D_p + (i * 2 * D_p);  % Centered at 1.5*D_p
        theta = linspace(-pi/2, pi/2, 50);  % Half-circle (from -pi/2 to pi/2)
        
        % Calculate the x-coordinate of the outer bend based on parallel lines
        x_bend = x_end + 1.5*D_p;  % Set the x-coordinate after the last parallel line
        
        % Tangent to x = 1 for even i
        x_outer = x_bend + r_outer * cos(theta);  % Use standard half-circle
        
        % Adjust x_outer based on the parity of i
        x_outer = x_outer - 1.5 * D_p;  % Add 1.5 * D_p
        
        y_outer = y_outer_center + r_outer * sin(theta);
        plot(x_outer, y_outer, 'k-', 'LineWidth', 2);
    end
    
    % Plot outer half-circle bends (centered at 1.5 * D_p) for odd indices
    for i = 1:2:N-1  % Odd indices only
        r_outer = 1.5 * D_p;
        y_outer_center = 1.5 * D_p + (i * 2 * D_p);  % Centered at 1.5*D_p
        theta = linspace(-pi/2, pi/2, 50);  % Half-circle (from -pi/2 to pi/2)
        
        % Calculate the x-coordinate of the outer bend based on parallel lines
        x_bend = x_end + 1.5*D_p;  % Set the x-coordinate after the last parallel line
        
        % Tangent to x = 0 for odd i
        x_outer = x_bend - r_outer * cos(theta);  % Mirror half-circle (reverse direction)
        
        % Adjust x_outer based on the parity of i
        x_outer = x_outer + 1.5 * D_p - L_e ;  % Subtract 1.5 * D_p
        
        y_outer = y_outer_center + r_outer * sin(theta);
        plot(x_outer, y_outer, 'k-', 'LineWidth', 2);
    end

    % Plot additional 2 segments starting at (0,0) and ending at (0, 1.5 * D_p) and (D_p, 1.5 * D_p)
    x_segment1 = [0, 1.5 * D_p]; % x-coordinates for the first segment
    y_segment1 = [0, 0]; % y-coordinates for the first segment
    plot(x_segment1, y_segment1, 'k-', 'LineWidth', 2);  % Plot the first segment

    x_segment2 = [0, 1.5 * D_p]; % x-coordinates for the second segment
    y_segment2 = [D_p, D_p]; % y-coordinates for the second segment
    plot(x_segment2, y_segment2, 'k-', 'LineWidth', 2);  % Plot the second segment

        
    if rem(N,2) == 0
        % Plot additional 2 segments starting at (0,0) and ending at (0, 1.5 * D_p) and (D_p, 1.5 * D_p)
        x_segment3 = [L_e-(1.5*D_p), L_e]; % x-coordinates for the first segment
        y_segment3 = [D_p * (2*N + 1), D_p * (2*N +1)]; % y-coordinates for the first segment
        plot(x_segment3, y_segment3, 'k-', 'LineWidth', 2);  % Plot the first segment
    
        x_segment4 = [L_e-(1.5*D_p), L_e]; % x-coordinates for the second segment
        y_segment4 = [D_p*(2*N),D_p*2*N]; % y-coordinates for the second segment
        plot(x_segment4, y_segment4, 'k-', 'LineWidth', 2);  % Plot the second segment
    else
        % Plot additional 2 segments starting at (0,0) and ending at (0, 1.5 * D_p) and (D_p, 1.5 * D_p)
        x_segment3 = [0, 1.5*D_p]; % x-coordinates for the first segment
        y_segment3 = [D_p * (2*N + 1), D_p * (2*N +1)]; % y-coordinates for the first segment
        plot(x_segment3, y_segment3, 'k-', 'LineWidth', 2);  % Plot the first segment
    
        x_segment4 = [0,1.5*D_p]; % x-coordinates for the second segment
        y_segment4 = [D_p*(2*N),D_p*(2*N)]; % y-coordinates for the second segment
        plot(x_segment4, y_segment4, 'k-', 'LineWidth', 2);  % Plot the second segment
    end

    hold off;
        
end
