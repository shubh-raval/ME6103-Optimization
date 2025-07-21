% Clear workspace and close figures
clear all;
close all;

% Load the test results using readtable instead of load
beam_test_results = readtable('beam_test_results.csv');




% Convert table to array for easier processing
results = table2array(beam_test_results);

% Build polynomial regression model from the test data
X = results(:, 2:6); % Features: a, b, w, t, load
y = results(:, 7);   % Target: deflection

gprMdl = fitrgp(X, y, 'KernelFunction', 'ardsquaredexponential', ...
                      'Standardize', true);
% Instead of using 'poly3' or 'poly5', use interactions with specified degree
% This creates a model with linear terms, squared terms, and interactions
poly_mdl = fitlm(X, y, 'purequadratic', 'RobustOpts', 'on');

% Define parameter bounds
h_bounds = [0.01, 0.2];        % height in meters
d_base_bounds = [0.095, 0.125]; % platform diameter in meters

% Create 2^2 factorial design with center point
design_points = [
    h_bounds(1), d_base_bounds(1); % low, low
    h_bounds(1), d_base_bounds(2); % low, high
    h_bounds(2), d_base_bounds(1); % high, low
    h_bounds(2), d_base_bounds(2); % high, high
    mean(h_bounds), mean(d_base_bounds) % center point
];

% Young's modulus in Pascals (aluminum)
E = 69e9;

% Calculate ellipse parameters for each design point
num_configs = size(design_points, 1);
beam_configs = zeros(num_configs, 4); % [a, b, w, t]

% Define fixed parameters for width and thickness
w_values = linspace(0.005, 0.05, num_configs); % width in meters
t_values = linspace(0.005, 0.05, num_configs); % thickness in meters

% Calculate ellipse parameters for each design point
for i = 1:num_configs
    h = design_points(i, 1);
    platform_diam = design_points(i, 2);
    
    % Calculate ellipse parameters
    b = platform_diam/2;
    base_rad = b; % Assuming base_rad = b = platform_diam/2
    a = sqrt(h^2 + base_rad^2);
    
    % Store parameters: a, b, w, t
    beam_configs(i, 1) = a;
    beam_configs(i, 2) = b;
    beam_configs(i, 3) = w_values(i);
    beam_configs(i, 4) = t_values(i);
end

% Define configuration names
config_names = cell(num_configs, 1);
for i = 1:num_configs
    if i < num_configs
        if design_points(i,1) == h_bounds(1)
            h_level = 'Low';
        else
            h_level = 'High';
        end
        
        if design_points(i,2) == d_base_bounds(1)
            d_level = 'Low';
        else
            d_level = 'High';
        end
        
        config_names{i} = sprintf('h: %s, d: %s', h_level, d_level);
    else
        config_names{i} = 'Center Point';
    end
end

% Define specific load points
loads = [100, 300, 500, 700];

% Process each beam configuration
for config_idx = 1:num_configs
    a = beam_configs(config_idx, 1);
    b = beam_configs(config_idx, 2);
    w = beam_configs(config_idx, 3);
    t = beam_configs(config_idx, 4);
    
    % Calculate second moment of area
    I = (w * t^3) / 12;
    
    % Initialize arrays for deflections
    meta_deflections = zeros(size(loads));
    cantilever_deflections = zeros(size(loads));
    poly_deflections = zeros(size(loads));
    
    % Calculate deflections for each load point
    for i = 1:length(loads)
        P = loads(i);
        
        % GPR meta-model prediction
        new_beam = [a, b, w, t, P];
        meta_deflections(i) = predict(gprMdl, new_beam);
        
        % Polynomial regression prediction
        poly_deflections(i) = predict(poly_mdl, new_beam);
        
        % Cantilever beam equation
        quarter_circumference = (pi/4) * (3*(a + b) - sqrt((3*a + b)*(a + 3*b)));
        cantilever_deflections(i) = P * quarter_circumference^3 / (3 * E * I) * 1000; % Convert to mm
    end
    
    % Plot the comparison for this configuration
    %plot_beam(loads, meta_deflections, cantilever_deflections, poly_deflections, config_names{config_idx});
end

% Helper function to plot beam deflection comparison
function plot_beam(loads, meta_deflections, cantilever_deflections, poly_deflections, config_name)
    figure;
    
    % Plot all three models
    scatter(loads, meta_deflections, 'b', 'filled');
    hold on;
    scatter(loads, cantilever_deflections, 'r', 'filled');
    scatter(loads, poly_deflections, 'g', 'filled');
    
    % Connect points with lines
    plot(loads, meta_deflections, 'b-');
    plot(loads, cantilever_deflections, 'r--');
    plot(loads, poly_deflections, 'g-.');
    
    % Add labels
    xlabel('Load (N)');
    ylabel('Deflection (mm)');
    title(['Beam Deflection: ', config_name]);
    legend('GPR Meta-model', 'Cantilever', 'Polynomial Regression');
    grid on;
end


% After the previous loop, add this section:

% Select 5 beams from the test data for comparison
num_test_beams = min(5, size(results, 1)); % In case there are fewer than 5 beams
test_beam_indices = round(linspace(1, size(results, 1), num_test_beams));

fprintf('Comparing models with %d actual test beams\n', num_test_beams);

% For each selected test beam
for beam_idx = 1:num_test_beams
    row_idx = test_beam_indices(beam_idx);
    
    % Extract beam parameters and actual deflection
    beam_params = results(row_idx, 2:6); % a, b, w, t, load
    actual_deflection = results(row_idx, 7);
    
    % Extract individual parameters
    a = beam_params(1);
    b = beam_params(2);
    w = beam_params(3);
    t = beam_params(4);
    P = beam_params(5); % Load
    
    % Calculate second moment of area
    I = (w * t^3) / 12;
    
    % Calculate deflections using all three models
    gpr_deflection = predict(gprMdl, beam_params);
    poly_deflection = predict(poly_mdl, beam_params);
    
    % Cantilever beam equation
    quarter_circumference = (pi/4) * (3*(a + b) - sqrt((3*a + b)*(a + 3*b)));
    cantilever_deflection = P * quarter_circumference^3 / (3 * E * I) * 1000; % Convert to mm
    
    % Create array of model predictions
    model_names = {'Actual', 'GPR Model', 'Poly Model', 'Cantilever'};
    deflections = [actual_deflection, gpr_deflection, poly_deflection, cantilever_deflection];
    
    % Calculate errors
    errors = deflections(2:end) - deflections(1);
    rel_errors = (errors / actual_deflection) * 100;
    
    fprintf('Beam %d (a=%.3f, b=%.3f, w=%.3f, t=%.3f, P=%.1f N):\n', beam_idx, a, b, w, t, P);
    fprintf('  Actual: %.7f mm\n', actual_deflection);
    fprintf('  GPR: %.7f mm (error: %.7f mm, %.7f%%)\n', gpr_deflection, errors(1), rel_errors(1));
    fprintf('  Poly: %.7f mm (error: %.7f mm, %.7f%%)\n', poly_deflection, errors(2), rel_errors(2));
    fprintf('  Cantilever: %.7f mm (error: %.7f mm, %.7f%%)\n', cantilever_deflection, errors(3), rel_errors(3));
    
    % Plot comparison for this test beam
    %plot_beam_comparison(model_names, deflections, beam_idx);
    %plot_all_model_predictions(model_names, deflections, beam_idx);

end

function plot_all_model_predictions(model_names, deflections, beam_idx)
    % Extract values
    actual = deflections(1);
    predictions = deflections(2:end); % GPR, Poly, Cantilever
    
    % Plot
    figure;
    hold on;

    % Plot each prediction
    colors = {'b', 'g', [0.8500 0.3250 0.0980]}; % blue, green, orange
    markers = {'o', 's', '^'};
    for i = 1:length(predictions)
        scatter(actual, predictions(i), 100, 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', colors{i}, ...
            'Marker', markers{i}, ...
            'DisplayName', model_names{i+1});
    end

    % Plot 1:1 reference line
    max_val = max([actual, predictions]) * 1.1;
    plot([0, max_val], [0, max_val], 'k--', 'LineWidth', 1.5, 'DisplayName', '1:1 Line');

    % Axis/label stuff
    xlabel('Actual Deflection (mm)');
    ylabel('Predicted Deflection (mm)');
    title(sprintf('Prediction Comparison for Test Beam %d', beam_idx));
    legend('Location', 'best');
    grid on;
    axis equal;
    xlim([0, max_val]);
    ylim([0, max_val]);

    hold off;
end


% Helper function to plot model comparison for test beams
function plot_beam_comparison(model_names, deflections, beam_idx)
    figure;
    
    % Create bar plot of deflections
    bar(deflections);
    
    % Set x-tick labels
    set(gca, 'XTickLabel', model_names);
    
    % Add labels
    ylabel('Deflection (mm)');
    title(sprintf('Model Comparison for Test Beam %d', beam_idx));
    
    % Add text labels with values above each bar
    for i = 1:length(deflections)
        text(i, deflections(i) + max(deflections)*0.03, sprintf('%.2f', deflections(i)), ...
             'HorizontalAlignment', 'center');
    end
    
    % Add grid
    grid on;
    
    % Optional: Add error percentage labels for model predictions
    for i = 2:length(deflections)
        error_pct = ((deflections(i) - deflections(1)) / deflections(1)) * 100;
        text(i, deflections(i) - max(deflections)*0.07, sprintf('%.1f%%', error_pct), ...
             'HorizontalAlignment', 'center', 'Color', 'red');
    end
end 

all_actual = [];
all_predicted = [];
all_models = [];
all_loads = [];
all_beam_ids = []; % NEW: Track beam numbers

% Use all 20 test beams from results
num_test_beams = 20;
test_beam_indices = 1:num_test_beams;

fprintf('Comparing models with %d actual test beams\n', num_test_beams);

for beam_idx = 1:num_test_beams
    row_idx = test_beam_indices(beam_idx);
    beam_params = results(row_idx, 2:6);
    actual_deflection = results(row_idx, 7);
    
    a = beam_params(1); b = beam_params(2);
    w = beam_params(3); t = beam_params(4); P = beam_params(5);
    I = (w * t^3) / 12;
    
    gpr_deflection = predict(gprMdl, beam_params);
    poly_deflection = predict(poly_mdl, beam_params);
    quarter_circ = (pi/4) * (3*(a + b) - sqrt((3*a + b)*(a + 3*b)));
    cantilever_deflection = P * quarter_circ^3 / (3 * E * I) * 1000;
    
    % Collect results
    all_actual     = [all_actual; repmat(actual_deflection, 3, 1)];
    all_predicted  = [all_predicted; gpr_deflection; poly_deflection; cantilever_deflection];
    all_models     = [all_models; "GPR"; "Poly"; "Cantilever"];
    all_loads      = [all_loads; repmat(P, 3, 1)];
    all_beam_ids   = [all_beam_ids; repmat(beam_idx, 3, 1)]; % NEW
end

% Call plotting function with beam IDs
plot_summary_by_load(all_actual, all_predicted, all_models, all_loads, all_beam_ids);


function plot_summary_by_load(all_actual, all_predicted, all_models, all_loads, all_beam_ids)
    unique_loads = unique(all_loads);
    models = unique(all_models);
    colors = lines(numel(models));
    markers = {'o', 's', '^'}; % Adjust as needed

    for i = 1:length(unique_loads)
        figure;
        load_val = unique_loads(i);
        idx_load = all_loads == load_val;

        for j = 1:length(models)
            model_name = models(j);
            model_idx = all_models == model_name & idx_load;

            scatter(all_actual(model_idx), all_predicted(model_idx), ...
                80, ...
                'Marker', markers{j}, ...
                'MarkerFaceColor', colors(j,:), ...
                'MarkerEdgeColor', 'k', ...
                'DisplayName', char(model_name));
            hold on;

            % Annotate with beam ID
            model_beam_ids = all_beam_ids(model_idx);
            actual_vals = all_actual(model_idx);
            predicted_vals = all_predicted(model_idx);
            for k = 1:length(model_beam_ids)
                text(actual_vals(k) + 0.02, predicted_vals(k), ...
                    sprintf('%d', model_beam_ids(k)), ...
                    'FontSize', 8, 'Color', 'k');
            end
        end

        % 1:1 reference line
        max_val = max([all_actual(idx_load); all_predicted(idx_load)]) * 1.1;
        plot([0, max_val], [0, max_val], 'k--', 'LineWidth', 1.2, 'DisplayName', '1:1 Line');

        title(sprintf('Prediction Comparison (Load = %.1f N)', load_val));
        xlabel('Actual Deflection (mm)');
        ylabel('Predicted Deflection (mm)');
        grid on;
        axis equal;
        xlim([0, max_val]);
        ylim([0, max_val]);
        legend('Location', 'best');
    end
end


function plot_summary_scatter(all_actual, all_predicted, all_models, all_loads)
    figure;
    hold on;

    % Unique models
    models = unique(all_models);
    colors = lines(numel(models)); % auto-colors
    markers = {'o', 's', '^'}; % You can extend this if more models

    for i = 1:numel(models)
        idx = all_models == models(i);
        scatter(all_actual(idx), all_predicted(idx), ...
            80, all_loads(idx), ... % Marker size/colors based on load
            'Marker', markers{i}, ...
            'MarkerFaceColor', colors(i,:), ...
            'MarkerEdgeColor', 'k', ...
            'DisplayName', models(i));
    end

    % 1:1 reference line
    max_val = max([all_actual; all_predicted]) * 1.1;
    plot([0, max_val], [0, max_val], 'k--', 'LineWidth', 1.2, 'DisplayName', '1:1 Line');

    xlabel('Actual Deflection (mm)');
    ylabel('Predicted Deflection (mm)');
    title('Model Prediction Summary for All Test Beams');
    colorbar;
    colormap(parula);
    c = colorbar;
    c.Label.String = 'Applied Load (N)';

    legend('Location', 'best');
    axis equal;
    grid on;
    xlim([0, max_val]);
    ylim([0, max_val]);

    hold off;
end
