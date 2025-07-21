% NSGA2_V2.m
% Enhanced NSGA-II implementation for SPM design optimization

function NSGA2_V2()
    %% Parameters
    params.load = 9.81;
    params.platform_diam = 0.10;
    params.base_rad = 0.16;
    params.h = 0.20;
    params.theta_deg = 15;
    params.actuator_masses = [41, 90, 107, 120, 122, 135, 195, 207.8, 230, 248, 248, 260, 282, 300, 450, 450, 570, 4000];
    params.debug = false;

    %% NSGA-II Enhanced Settings
    pop_size = 50;       % Increased population size
    max_gens = 30;       % More generations for exploration
    num_vars = 3;        % Number of decision variables

    % Decision variable bounds
    lb = [0.005, 0.005, 1]; 
    ub = [0.05, 0.05, length(params.actuator_masses)];

    % Advanced NSGA-II parameters
    crossover_prob = 0.9;        % Higher crossover probability
    mutation_prob = 2 / num_vars;% Increased mutation probability
    eta_c = 30;                  % Higher distribution index for crossover
    eta_m = 30;                  % Higher distribution index for mutation

    disp('Starting Enhanced NSGA-II optimization');

    %% Initialize Population with Advanced Diversity Strategy
    pop = zeros(pop_size, num_vars);

    % Create diverse initial population with advanced stratification
    % Group 1: Lightweight designs (33%)
    for i = 1:floor(pop_size/3)
        pop(i,1) = lb(1) + (ub(1)-lb(1)) * rand(); % Full range for W
        pop(i,2) = lb(2) + (ub(2)-lb(2)) * rand(); % Full range for T
        pop(i,3) = round(1 + (length(params.actuator_masses)/3)*rand()); % Lighter actuators
    end

    % Group 2: Mid-range designs (33%)
    for i = floor(pop_size/3)+1:2*floor(pop_size/3)
        pop(i,1) = lb(1) + (ub(1)-lb(1)) * 0.5 + (ub(1)-lb(1)) * 0.25 * rand(); % Centered with more spread
        pop(i,2) = lb(2) + (ub(2)-lb(2)) * 0.5 + (ub(2)-lb(2)) * 0.25 * rand(); % Centered with more spread
        pop(i,3) = round((length(params.actuator_masses)/3) + (length(params.actuator_masses)/3)*rand()); % Mid-range actuators
    end

    % Group 3: Exploratory designs (33%)
    for i = 2*floor(pop_size/3)+1:pop_size
        pop(i,1) = lb(1) + (ub(1)-lb(1)) * rand(); % Full range exploration
        pop(i,2) = lb(2) + (ub(2)-lb(2)) * rand(); % Full range exploration
        pop(i,3) = round(1 + (length(params.actuator_masses)-1)*rand()); % Any actuator
    end

    % Ensure all actuator indices are valid
    for i = 1:pop_size
        pop(i,3) = max(1, min(length(params.actuator_masses), round(pop(i,3))));
    end

    % Evaluate initial population
    F = zeros(pop_size, 2);
    disp('Evaluating initial population...');
    for i = 1:pop_size
        W = pop(i,1);
        T = pop(i,2);
        actuator_idx = round(pop(i,3));
        x = [W, T, actuator_idx];
        F(i,:) = SPM_objectives(x, params);
        
        % Print out to monitor diversity
        valid_str = '-valid-';
        if F(i,2) > -0.5
            valid_str = '-invalid-';
        end
        
        disp(['Individual ', num2str(i), ': W=', num2str(W), ', T=', num2str(T), ...
              ', Act=', num2str(actuator_idx), ' -> Mass=', num2str(F(i,1)), ...
              'g, GCI=', num2str(-F(i,2)), ' ', valid_str]);
    end

    %% Main NSGA-II Loop with Enhanced Exploration
    best_mass = Inf;
    best_gci = -Inf;
    best_design = zeros(1, num_vars);

    % Store history data
    history.gen = 1:max_gens;
    history.best_mass = zeros(max_gens, 1);
    history.best_gci = zeros(max_gens, 1);
    history.avg_mass = zeros(max_gens, 1);
    history.avg_gci = zeros(max_gens, 1);
    history.pop = cell(max_gens, 1);
    history.F = cell(max_gens, 1);

    % Custom NSGA-II with advanced diversity preservation
    for gen = 1:max_gens
        disp(['Generation ', num2str(gen), '/', num2str(max_gens)]);
        
        % Generate Offspring
        offspring = zeros(pop_size, num_vars);
        for i = 1:2:pop_size
            % Select two parents
            p1 = tournament_selection(F, pop_size);
            p2 = tournament_selection(F, pop_size);
            
            % Apply crossover
            if rand < crossover_prob
                [child1, child2] = sbx_crossover(pop(p1,:), pop(p2,:), lb, ub, eta_c);
            else
                child1 = pop(p1,:);
                child2 = pop(p2,:);
            end
            
            % Apply mutation
            child1 = polynomial_mutation(child1, lb, ub, mutation_prob, eta_m);
            child2 = polynomial_mutation(child2, lb, ub, mutation_prob, eta_m);
            
            % Discretize actuator selection
            child1(3) = round(child1(3));
            child2(3) = round(child2(3));
            
            % Ensure actuator indices are valid
            child1(3) = max(1, min(length(params.actuator_masses), child1(3)));
            child2(3) = max(1, min(length(params.actuator_masses), child2(3)));
            
            % Add to offspring
            if i < pop_size
                offspring(i,:) = child1;
                offspring(i+1,:) = child2;
            else
                offspring(i,:) = child1;
            end
        end
        
        % Evaluate offspring
        F_off = zeros(pop_size, 2);
        disp('Evaluating offspring...');
        for i = 1:pop_size
            W = offspring(i,1);
            T = offspring(i,2);
            actuator_idx = round(offspring(i,3));
            
            x = [W, T, actuator_idx];
            F_off(i,:) = SPM_objectives(x, params);
        end
        
        % Combine parent and offspring populations
        combined_pop = [pop; offspring];
        combined_F = [F; F_off];
        
        % Non-dominated sorting
        [fronts, ~] = non_dominated_sort(combined_F);
        
        % Calculate crowding distance
        distances = crowding_distance(fronts, combined_F);
        
        % Select next generation with diversity preservation
        new_pop = [];
        new_F = [];
        front_idx = 1;
        
        % Select from each front, preserving diversity
        while size(new_pop, 1) + length(fronts{front_idx}) <= pop_size
            front = fronts{front_idx};
            new_pop = [new_pop; combined_pop(front, :)];
            new_F = [new_F; combined_F(front, :)];
            front_idx = front_idx + 1;
            
            if front_idx > length(fronts)
                break;
            end
        end
        
        % Handle the last front - try to maintain diversity
        if size(new_pop, 1) < pop_size && front_idx <= length(fronts)
            last_front = fronts{front_idx};
            
            % Sort primarily by crowding distance
            [~, sorted_idx] = sort(distances(last_front), 'descend');
            last_front_sorted = last_front(sorted_idx);
            
            % Calculate how many more individuals we need
            remaining = pop_size - size(new_pop, 1);
            
            % Select remaining individuals, but prioritize diversity
            if remaining > 0
                % Add some from the last front
                selected = last_front_sorted(1:min(remaining, length(last_front_sorted)));
                new_pop = [new_pop; combined_pop(selected, :)];
                new_F = [new_F; combined_F(selected, :)];
            end
        end
        
        % If we still don't have enough, add random points
        while size(new_pop, 1) < pop_size
            % Generate a random individual
            random_ind = lb + rand(1, num_vars) .* (ub - lb);
            random_ind(3) = round(random_ind(3));
            random_ind(3) = max(1, min(length(params.actuator_masses), random_ind(3)));
            
            % Evaluate
            x = [random_ind(1), random_ind(2), round(random_ind(3))];
            random_f = SPM_objectives(x, params);
            
            % Add to population
            new_pop = [new_pop; random_ind];
            new_F = [new_F; random_f];
        end
        
        % Periodically inject diversity (every 5 generations)
        if mod(gen, 5) == 0
            % Inject some random solutions to avoid premature convergence
            num_random = floor(pop_size * 0.1);
            
            if num_random > 0
                % Generate random solutions across the range
                for i = 1:num_random
                    rand_idx = randi(pop_size);
                    
                    % Create a random solution
                    if rand < 0.5
                        % Lightweight solution
                        new_pop(rand_idx, 1) = lb(1) + (ub(1)-lb(1))/3 * rand();
                        new_pop(rand_idx, 2) = lb(2) + (ub(2)-lb(2))/3 * rand();
                        new_pop(rand_idx, 3) = round(1 + 5*rand()); % Lighter actuators
                    else
                        % Heavier solution
                        new_pop(rand_idx, 1) = lb(1) + (ub(1)-lb(1))/2 + (ub(1)-lb(1))/2 * rand();
                        new_pop(rand_idx, 2) = lb(2) + (ub(2)-lb(2))/2 + (ub(2)-lb(2))/2 * rand();
                        new_pop(rand_idx, 3) = round(1 + (length(params.actuator_masses)-1)*rand());
                    end
                    
                    % Validate actuator index
                    new_pop(rand_idx, 3) = max(1, min(length(params.actuator_masses), round(new_pop(rand_idx, 3))));
                    
                    % Re-evaluate this individual
                    x = [new_pop(rand_idx, 1), new_pop(rand_idx, 2), new_pop(rand_idx, 3)];
                    new_F(rand_idx, :) = SPM_objectives(x, params);
                end
                
                disp(['Injected ', num2str(num_random), ' random solutions for diversity']);
            end
        end
        
        % Update population
        pop = new_pop;
        F = new_F;
        
        % Update best solutions
        [min_mass, mass_idx] = min(F(:,1));
        [min_neg_gci, gci_idx] = min(F(:,2));
        max_gci = -min_neg_gci;
        
        if min_mass < best_mass
            best_mass = min_mass;
        end
        
        if max_gci > best_gci
            best_gci = max_gci;
            best_design = pop(gci_idx,:);
        end
        
        % Store history
        history.best_mass(gen) = best_mass;
        history.best_gci(gen) = best_gci;
        history.avg_mass(gen) = mean(F(:,1));
        history.avg_gci(gen) = mean(-F(:,2));
        history.pop{gen} = pop;
        history.F{gen} = F;
        
        % Display progress
        disp(['Best mass: ', num2str(best_mass), ' g, Best GCI: ', num2str(best_gci)]);
        
        % Check diversity metrics
        valid_designs = sum(-F(:,2) > 0.5);
        mass_range = max(F(:,1)) - min(F(:,1));
        disp(['  Valid designs: ', num2str(valid_designs), '/', num2str(pop_size), ...
              ', Mass range: ', num2str(mass_range), ' g']);
        
        % Save intermediate results
        if mod(gen, 5) == 0
            save(['nsga2_gen_', num2str(gen), '.mat'], 'pop', 'F', 'best_mass', 'best_gci', 'best_design', 'history');
        end
    end

    %% Final Results
    disp(' ');
    disp('=====================');
    disp('Optimization complete');
    disp('=====================');
    disp(' ');

    % Find valid designs
    valid_idx = find(-F(:,2) > 0.5);
    invalid_idx = find(-F(:,2) <= 0.5);

    disp(['Found ', num2str(length(valid_idx)), ' valid designs and ', ...
          num2str(length(invalid_idx)), ' invalid designs.']);

    % Top valid solutions by GCI
    disp('Top valid solutions by GCI:');
    [~, idx] = sort(F(valid_idx,2));
    for i = 1:min(5, length(idx))
        design_idx = valid_idx(idx(i));
        disp(['#', num2str(i), ': W=', num2str(pop(design_idx,1), '%.4f'), ...
              ', T=', num2str(pop(design_idx,2), '%.4f'), ...
              ', Actuator=', num2str(round(pop(design_idx,3))), ...
              ' -> Mass=', num2str(F(design_idx,1), '%.2f'), ...
              ' g, GCI=', num2str(-F(design_idx,2), '%.6f')]);
    end

    % Top valid solutions by Mass
    disp(' ');
    disp('Top valid solutions by Mass:');
    [~, idx] = sort(F(valid_idx,1));
    for i = 1:min(5, length(idx))
        design_idx = valid_idx(idx(i));
        disp(['#', num2str(i), ': W=', num2str(pop(design_idx,1), '%.4f'), ...
              ', T=', num2str(pop(design_idx,2), '%.4f'), ...
              ', Actuator=', num2str(round(pop(design_idx,3))), ...
              ' -> Mass=', num2str(F(design_idx,1), '%.2f'), ...
              ' g, GCI=', num2str(-F(design_idx,2), '%.6f')]);
    end

    % Best design information
    disp(' ');
    disp('Best design (highest GCI):');
    disp(['W=', num2str(best_design(1), '%.4f'), ...
          ', T=', num2str(best_design(2), '%.4f'), ...
          ', Actuator=', num2str(round(best_design(3)))]);
    disp(['Mass=', num2str(best_mass, '%.2f'), ' g, GCI=', num2str(best_gci, '%.6f')]);

    %% Visualization
    
    % 1. PARETO FRONT WITH VALID/INVALID CLASSIFICATION
    figure('Position', [100, 100, 900, 700]);
    hold on;

    % Plot invalid designs
    if ~isempty(invalid_idx)
        scatter(F(invalid_idx,1), -F(invalid_idx,2), 70, 'ro', 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', 'r');
    end

    % Plot valid designs
    if ~isempty(valid_idx)
        scatter(F(valid_idx,1), -F(valid_idx,2), 100, 'bo', 'filled');
        
        % Connect valid designs to show the Pareto front
        [sorted_mass_valid, sort_idx] = sort(F(valid_idx,1));
        plot(sorted_mass_valid, -F(valid_idx(sort_idx),2), 'b-', 'LineWidth', 2);
    end

    % Highlight best design
    [~, best_idx] = max(-F(:,2));
    scatter(F(best_idx,1), -F(best_idx,2), 200, 'gd', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    text(F(best_idx,1)+5, -F(best_idx,2), 'Best Design', 'FontSize', 12, 'FontWeight', 'bold');

    % Add design parameters for selected points
    for i = 1:length(valid_idx)
        if mod(i, 3) == 0 || i == 1 || i == length(valid_idx)
            idx = valid_idx(i);
            text(F(idx,1)+2, -F(idx,2), sprintf('W=%.3f\nT=%.3f\nA=%d', ...
                 pop(idx,1), pop(idx,2), round(pop(idx,3))), 'FontSize', 8);
        end
    end

    % Set labels and grid
    grid on;
    xlabel('Total Mass [g]', 'FontSize', 14);
    ylabel('Global Conditioning Index (GCI)', 'FontSize', 14);
    title('Pareto Front with Valid/Invalid Designs', 'FontSize', 16);
    legend('Invalid Designs', 'Valid Designs', 'Pareto Front', 'Best Design', 'Location', 'best');

    % 2. MASS VS DESIGN VARIABLES COLORED BY GCI
    figure('Position', [100, 100, 1000, 500]);

    % Plot W vs Mass colored by GCI
    subplot(1,2,1);
    scatter(pop(:,1), F(:,1), 80, -F(:,2), 'filled');
    hold on;
    scatter(pop(valid_idx,1), F(valid_idx,1), 100, -F(valid_idx,2), 'filled', 'MarkerEdgeColor', 'k');
    colormap(jet);
    colorbar;
    xlabel('Link Width (W) [m]', 'FontSize', 12);
    ylabel('Mass [g]', 'FontSize', 12);
    title('Mass vs. Link Width', 'FontSize', 14);
    grid on;

    % Plot T vs Mass colored by GCI
    subplot(1,2,2);
    scatter(pop(:,2), F(:,1), 80, -F(:,2), 'filled');
    hold on;
    scatter(pop(valid_idx,2), F(valid_idx,1), 100, -F(valid_idx,2), 'filled', 'MarkerEdgeColor', 'k');
    colormap(jet);
    c = colorbar;
    ylabel(c, 'GCI Value', 'FontSize', 12);
    xlabel('Link Thickness (T) [m]', 'FontSize', 12);
    ylabel('Mass [g]', 'FontSize', 12);
    title('Mass vs. Link Thickness', 'FontSize', 14);
    grid on;

    % 3. DESIGN ACTUATOR ANALYSIS
    figure('Position', [100, 100, 1200, 500]);

    % Get unique actuator types
    unique_actuators = unique(round(pop(:,3)));

    % Subplot for Mass distribution
    subplot(1,2,1);
    boxplot(F(:,1), round(pop(:,3)), 'Labels', arrayfun(@num2str, unique_actuators, 'UniformOutput', false));
    title('Mass Distribution by Actuator Type', 'FontSize', 14);
    xlabel('Actuator Index', 'FontSize', 12);
    ylabel('Mass [g]', 'FontSize', 12);
    grid on;

    % Subplot for GCI distribution
    subplot(1,2,2);
    boxplot(-F(:,2), round(pop(:,3)), 'Labels', arrayfun(@num2str, unique_actuators, 'UniformOutput', false));
    title('GCI Distribution by Actuator Type', 'FontSize', 14);
    xlabel('Actuator Index', 'FontSize', 12);
    ylabel('GCI', 'FontSize', 12);
    grid on;

    % Save final results
    save('nsga2_final_results.mat', 'pop', 'F', 'best_mass', 'best_gci', 'best_design', 'history', 'valid_idx', 'invalid_idx');
end

%% Helper Functions

function idx = tournament_selection(F, pop_size)
    % Binary tournament selection
    i1 = randi(pop_size);
    i2 = randi(pop_size);
    
    % Check dominance
    if dominates(F(i1,:), F(i2,:))
        idx = i1;
    elseif dominates(F(i2,:), F(i1,:))
        idx = i2;
    else
        % No dominance, choose randomly
        if rand < 0.5
            idx = i1;
        else
            idx = i2;
        end
    end
end

function flag = dominates(x, y)
    % Check if x dominates y for minimization
    flag = all(x <= y) && any(x < y);
end

function [child1, child2] = sbx_crossover(parent1, parent2, lb, ub, eta)
    % Simulated Binary Crossover
    child1 = parent1;
    child2 = parent2;
    
    for i = 1:length(parent1)
        if abs(parent1(i) - parent2(i)) > 1e-10
            if parent1(i) < parent2(i)
                y1 = parent1(i);
                y2 = parent2(i);
            else
                y1 = parent2(i);
                y2 = parent1(i);
            end
            
            % Calculate offspring values
            r = rand();
            beta = 1.0 + (2.0 * (y1 - lb(i)) / (y2 - y1));
            alpha = 2.0 - beta^(-(eta + 1.0));
            
            if r <= 1.0 / alpha
                beta_q = (r * alpha)^(1.0 / (eta + 1.0));
            else
                beta_q = (1.0 / (2.0 - r * alpha))^(1.0 / (eta + 1.0));
            end
            
            c1 = 0.5 * ((y1 + y2) - beta_q * (y2 - y1));
            c2 = 0.5 * ((y1 + y2) + beta_q * (y2 - y1));
            
            c1 = max(lb(i), min(ub(i), c1));
            c2 = max(lb(i), min(ub(i), c2));
            
            if parent1(i) < parent2(i)
                child1(i) = c1;
                child2(i) = c2;
            else
                child1(i) = c2;
                child2(i) = c1;
            end
        end
    end
end

function child = polynomial_mutation(parent, lb, ub, prob, eta)
    % Polynomial mutation
    child = parent;
    
    for i = 1:length(parent)
        if rand < prob
            y = parent(i);
            delta_max = ub(i) - lb(i);
            
            r = rand();
            if r < 0.5
                delta = (2.0 * r)^(1.0 / (eta + 1.0)) - 1.0;
            else
                delta = 1.0 - (2.0 * (1.0 - r))^(1.0 / (eta + 1.0));
            end
            
            child(i) = y + delta * delta_max;
            child(i) = max(lb(i), min(ub(i), child(i)));
        end
    end
end

function [fronts, ranks] = non_dominated_sort(F)
    n = size(F, 1);
    fronts = {};
    ranks = zeros(n, 1);
    
    % Calculate domination for each solution
    S = cell(n, 1);     % Set of solutions dominated by i
    n_dom = zeros(n, 1); % Number of solutions dominating i
    
    for i = 1:n
        S{i} = [];
        
        for j = 1:n
            if i ~= j
                if dominates(F(i,:), F(j,:))
                    S{i} = [S{i} j];
                elseif dominates(F(j,:), F(i,:))
                    n_dom(i) = n_dom(i) + 1;
                end
            end
        end
        
        % First front
        if n_dom(i) == 0
            ranks(i) = 1;
            if isempty(fronts)
                fronts{1} = i;
            else
                fronts{1} = [fronts{1} i];
            end
        end
    end
    
    % Find subsequent fronts
    k = 1;
    while k <= length(fronts)
        Q = [];
        for i = fronts{k}
            for j = S{i}
                n_dom(j) = n_dom(j) - 1;
                if n_dom(j) == 0
                    ranks(j) = k + 1;
                    Q = [Q j];
                end
            end
        end
        
        if ~isempty(Q)
            fronts{k+1} = Q;
        end
        
        k = k + 1;
    end
end

function distances = crowding_distance(fronts, F)
    n = size(F, 1);
    distances = zeros(n, 1);
    
    % Process each front
    for k = 1:length(fronts)
        front = fronts{k};
        if length(front) <= 2
            distances(front) = Inf;
            continue;
        end
        
        % Process each objective
        for m = 1:size(F, 2)
            % Sort by current objective
            [sorted_obj, sorted_idx] = sort(F(front, m));
            sorted_idx = front(sorted_idx);
            
            % Set extremes to infinity
            distances(sorted_idx(1)) = Inf;
            distances(sorted_idx(end)) = Inf;
            
            % Calculate crowding distance
            for i = 2:length(sorted_idx)-1
                if (sorted_obj(end) - sorted_obj(1)) > 1e-10
                    distances(sorted_idx(i)) = distances(sorted_idx(i)) + ...
                        (sorted_obj(i+1) - sorted_obj(i-1)) / (sorted_obj(end) - sorted_obj(1));
                end
            end
        end
    end
end
    