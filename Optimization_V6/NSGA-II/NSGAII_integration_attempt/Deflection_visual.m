function beam_deflection_visualization()
    % Beam Deflection Geometric Visualization
    
    % Fixed parameters
    h = 0.2;  % Total height of the beam
    
    % Angles to visualize
    angles = [0, 15, 30, 45, 60];
    
    % Create a figure
    figure('Position', [100, 100, 1200, 800]);
    
    % Subplot for ellipse visualization
    subplot(2,2,[1,2]);
    hold on;
    
    % Color map for different angles
    colors = jet(length(angles));
    
    % Plot ellipses for different angles
    for i = 1:length(angles)
        theta = angles(i);
        
        % Calculate semi-major and semi-minor axes
        a = h;  % Major axis always full height
        b = h * cosd(theta);  % Minor axis reduces with angle
        
        % Create points for ellipse
        t = linspace(0, 2*pi, 100);
        x = a * cos(t);
        y = b * sin(t);
        
        % Plot the ellipse
        plot(x, y, 'Color', colors(i,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('θ = %d°', theta));
    end
    
    title('Beam Deflection: Ellipse Shape Variation', 'FontSize', 14);
    xlabel('Width', 'FontSize', 12);
    ylabel('Height', 'FontSize', 12);
    axis equal;
    grid on;
    legend('show', 'Location', 'best');
    
    % Subplot for axis length comparison
    subplot(2,2,3);
    major_lengths = ones(size(angles)) * h;
    minor_lengths = h * cosd(angles);
    
    bar([major_lengths; minor_lengths]');
    title('Axis Lengths vs. Angle', 'FontSize', 14);
    xlabel('Angle (degrees)', 'FontSize', 12);
    ylabel('Axis Length', 'FontSize', 12);
    legend('Major Axis', 'Minor Axis', 'Location', 'best');
    set(gca, 'XTickLabel', arrayfun(@num2str, angles, 'UniformOutput', false));
    
    % Subplot for area calculation
    subplot(2,2,4);
    areas = pi * h * (h * cosd(angles));
    
    plot(angles, areas, 'ro-', 'LineWidth', 2);
    title('Ellipse Area vs. Angle', 'FontSize', 14);
    xlabel('Angle (degrees)', 'FontSize', 12);
    ylabel('Ellipse Area', 'FontSize', 12);
    grid on;
    
    % Suptitle
    sgtitle('Geometric Interpretation of Beam Deflection', 'FontSize', 16);
    
    % Print out some numerical values for reference
    fprintf('Numerical Analysis of Beam Deflection Geometry:\n');
    fprintf('Fixed Height (h): %.4f m\n', h);
    fprintf('\nAngle | Major Axis | Minor Axis | Area\n');
    fprintf('----------------------------------\n');
    for i = 1:length(angles)
        theta = angles(i);
        major = h;
        minor = h * cosd(theta);
        area = pi * major * minor;
        fprintf('%5d° | %10.4f | %10.4f | %8.4f\n', theta, major, minor, area);
    end
end