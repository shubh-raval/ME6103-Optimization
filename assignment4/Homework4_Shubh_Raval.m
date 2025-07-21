%% Homework 4 Shubh Raval


% Create points
x1_range = linspace(0, 6, 100);
x2_range = linspace(0, 6, 100);
[X1, X2] = meshgrid(x1_range, x2_range);

% Create the Function 
F = @(x1,x2) cos(x1) + sin(x2);
Fplot = cos(X1) + sin(X2);

% Create 3D surface plot of the function 
figure;
surf(X1, X2, Fplot);
colormap('jet');
colorbar;
title('f(x_1, x_2) = cos(x_1) + sin(x_2)');
xlabel('x_1');
ylabel('x_2');
zlabel('f(x_1, x_2)');
grid on;
view(45, 30);
create = true;
true_func = @(x) F(x(1), x(2));
fixed_x2 = 3;  % for error checking
%% Create Uniform Grid of 25 Test Points 
lb = [0, 0]; 
ub = [6, 6]; 
[X1_grid, X2_grid] = meshgrid(linspace(lb(1), ub(1), 5), linspace(lb(2), ub(2), 5));
test_points = [X1_grid(:), X2_grid(:)];
test_responses = zeros(25, 1);
for i = 1:25
    test_responses(i) = F(test_points(i,1), test_points(i,2));
end

%% PROBLEM 1
disp('CREATE EXPERIMENTAL DATA SETS FOR PROBLEM 1')

[ccd, f_eval_ccd,halton20, f_eval_halton20] = createExperimentalsets(create);


%% PROBLEM 2
disp('BEGIN META MODELING FOR PROBLEM 2')


%% a. Second Order Regression Model based on CCD Data 
mdl = fitlm(ccd, f_eval_ccd, 'quadratic');

% Create prediction function from the model
predict_fitlm = @(x) predict(mdl, x);

[X1_vis, X2_vis] = meshgrid(linspace(lb(1), ub(1), 50), linspace(lb(2), ub(2), 50));
reg_fitlm_vis = zeros(50, 50);
for i = 1:50
    for j = 1:50
        reg_fitlm_vis(i,j) = predict_fitlm([X1_vis(i,j), X2_vis(i,j)]);
    end
end

ccd_predictions_fitlm = predict(mdl, ccd);

plotDesignPoints(X1_vis, X2_vis, reg_fitlm_vis, ccd, ccd_predictions_fitlm, 'a. CCD Second Order Regression');
[RAE_ccd_fitlm, RME_ccd_fitlm] = evaluateMetamodel(predict_fitlm, test_points, test_responses, fixed_x2, true_func);

%% b. Second Order Regression Model based on Halton-20 Data using fitlm
h_20_mdl = fitlm(halton20, f_eval_halton20, "quadratic");

% Create prediction function from the model
predict_h_20 = @(x) predict(h_20_mdl, x);

reg_fith20_vis = zeros(50, 50);
for i = 1:50
    for j = 1:50
        reg_fith20_vis(i,j) = predict_h_20([X1_vis(i,j), X2_vis(i,j)]);
    end
end

h_20_predictions = predict(h_20_mdl, halton20);

plotDesignPoints(X1_vis, X2_vis, reg_fith20_vis, halton20, h_20_predictions, 'b. Second Order Regression Model (Halton-20 Data)');

[RAE_h20_fitlm, RME_h20_fitlm] = evaluateMetamodel(predict_h_20, test_points, test_responses, fixed_x2,true_func);

%% c. Kriging Model based on CCD Data using 0 order polynomial and Gaussian Fit
d_ccdmodel = buildKrigingModel(ccd, f_eval_ccd, X1_vis, X2_vis);
predict_kriging_ccd = @(x) predictor(x, d_ccdmodel);
[RAE_ccd_kriging, RME_ccd_kriging] = evaluateMetamodel(predict_kriging_ccd, test_points, test_responses,fixed_x2,true_func);

%% d. Kriging Model based on Halton 20 Data using 0 order polynomial and Gaussian Fit
d_h_20model = buildKrigingModel(halton20, f_eval_halton20, X1_vis, X2_vis);
predict_kriging_h20 = @(x) predictor(x, d_h_20model);
[RAE_h20_kriging, RME_h20_kriging] = evaluateMetamodel(predict_kriging_h20, test_points, test_responses,fixed_x2,true_func);


function [RAAE, RMAE] = evaluateMetamodel(predict_func, test_points, test_responses, fixed_x2,true_func)

    test_predictions = zeros(size(test_points,1),1);
    for i = 1:size(test_points,1)
        test_predictions(i) = predict_func(test_points(i,:));
    end
    
    % Calculate errors
    errors = test_predictions - test_responses;
    abs_errors = abs(errors);
    
    response_std = std(test_responses);
    rel_errors = abs_errors / response_std;

    RAAE = mean(rel_errors);  % Relative Average Error
    RMAE = max(rel_errors);   % Relative Maximum Error
    

    fprintf('\n======= Metamodel Accuracy Metrics =======\n');
    disp(RAAE);
    disp(RMAE);


    x1_range = linspace(min(test_points(:,1)), max(test_points(:,1)), 200)';
    x2_const = fixed_x2 * ones(size(x1_range));
    input_grid = [x1_range, x2_const];
    true_vals = arrayfun(@(i) true_func(input_grid(i,:)), 1:length(x1_range))';
    pred_vals = arrayfun(@(i) predict_func(input_grid(i,:)), 1:length(x1_range))';

    % Plotting
    figure;
    plot(x1_range, true_vals, 'b-', 'LineWidth', 2); hold on;
    plot(x1_range, pred_vals, 'r--', 'LineWidth', 2);
    legend('True Function', 'Metamodel Prediction', 'Location', 'Best');
    xlabel('x1');
    ylabel('Function Value');
    title(sprintf('Metamodel vs True Function (x2 = %.2f)', fixed_x2));
    grid on;
    
end
function plotDesignPoints(X1, X2, Fplot, designPoints, f_eval, title_text)
    figure;
    
    % First subplot: 3D view
    subplot(1, 2, 1);
    surf(X1, X2, Fplot, 'FaceAlpha', 0.5);
    hold on;
    scatter3(designPoints(:,1), designPoints(:,2), f_eval, 100, 'r', 'filled');
    colormap('jet');
    title([title_text, ' - 3D View']);
    xlabel('x_1');
    ylabel('x_2');
    zlabel('f(x_1, x_2)');
    grid on;
    view(45, 30);
    
    % Second subplot: 2D view
    subplot(1, 2, 2);
    scatter(designPoints(:,1), designPoints(:,2), 80, f_eval, 'filled');
    colorbar;
    title([title_text, ' - 2D View']);
    xlabel('x_1');
    ylabel('x_2');
    axis([0 6 0 6]);
    grid on;
    
end
function natural = scaleFromCoded(coded, lb, ub)
    natural = lb + (coded + 1) .* (ub - lb) / 2;
end
function [ccd, f_eval_ccd,halton20, f_eval_halton20] = createExperimentalsets(create)
    if create
        % Define the domain boundaries
        x1_range = linspace(0, 6, 100);
        x2_range = linspace(0, 6, 100);
        
        % Create meshgrid for 3D surface plot
        [X1, X2] = meshgrid(x1_range, x2_range);
        
        % Create the Function to be Called in the metamodeling
        F = @(x1,x2) cos(x1) + sin(x2);
        
        % Compute function values
        Fplot = cos(X1) + sin(X2);
        %% a. Two-Level Factorial Design
        lb = [0, 0]; 
        ub = [6, 6];         
        factorialPoints = [lb(1), lb(2);
                          lb(1), ub(2);
                          ub(1), lb(2);
                          ub(1), ub(2)];
        
        f_eval_factorial = zeros(size(factorialPoints,1),1);
        for i = 1:size(factorialPoints,1)
            f_eval_factorial(i) = F(factorialPoints(i,1), factorialPoints(i,2));
        end
        
        plotDesignPoints(X1, X2, Fplot, factorialPoints, f_eval_factorial, 'a. Two-Level Factorial Design');
        
        %% b. Central Composite Design (CCD)
        alpha = sqrt(2);  
        axial_coded = [ alpha,  0;
                       -alpha,  0;
                        0,    alpha;
                        0,   -alpha];
        ccd_axial = scaleFromCoded(axial_coded, lb, ub);
        ccd_center = (lb + ub) / 2;
        ccd = [factorialPoints; ccd_axial; ccd_center];
        f_eval_ccd = zeros(size(ccd,1), 1);
        for i = 1:size(ccd, 1)
            f_eval_ccd(i) = F(ccd(i,1), ccd(i,2));
        end
        
        plotDesignPoints(X1, X2, Fplot, ccd, f_eval_ccd, 'CCD');
        
        %% c. Latin Hypercube Design with 10 Points
        lhsPoints = lhsdesign(10, 2);  
        % Need to scale the points 
        lhsPoints = lb + lhsPoints .* (ub - lb);
        
        f_eval_lhs = zeros(size(lhsPoints,1),1);
        for i = 1:size(lhsPoints,1)
            f_eval_lhs(i) = F(lhsPoints(i,1), lhsPoints(i,2));
        end
        
        % Plot the Latin Hypercube design
        plotDesignPoints(X1, X2, Fplot, lhsPoints, f_eval_lhs, 'c. Latin Hypercube Design (10 Points)');
        
        %% d. Halton Sequence Designs
        % For 10 points:
        pHalton10 = haltonset(2);
        halton10 = net(pHalton10,10);
        % Scale like lhs
        halton10 = lb + halton10 .* (ub - lb);
        
        f_eval_halton10 = zeros(size(halton10,1),1);
        for i = 1:size(halton10,1)
            f_eval_halton10(i) = F(halton10(i,1), halton10(i,2));
        end
        
        % Plot the Halton-10 design
        plotDesignPoints(X1, X2, Fplot, halton10, f_eval_halton10, 'd. Halton Sequence Design (10 Points)');
        
        % For 20 points:
        pHalton20 = haltonset(2);
        halton20 = net(pHalton20,20);

        % Scale like lhs
        halton20 = lb + halton20 .* (ub - lb);
        
        f_eval_halton20 = zeros(size(halton20,1),1);
        for i = 1:size(halton20,1)
            f_eval_halton20(i) = F(halton20(i,1), halton20(i,2));
        end
        
        % Plot the Halton-20 design
        plotDesignPoints(X1, X2, Fplot, halton20, f_eval_halton20, 'd. Halton Sequence Design (20 Points)');
    end
end
function dmodel = buildKrigingModel(design_points, f_eval, X1_vis, X2_vis)

    % Create the model
    xTrain = design_points;       
    yTrain = f_eval; 

    [xUnique, ia, ~] = unique(xTrain, 'rows', 'stable');
    yUnique = yTrain(ia);
    regr   = @regpoly0;   
    corr   = @corrgauss;  
    theta0 = 0.1;
    lob    = 1e-6;
    upb    = 10;

    [dmodel, ~] = dacefit(xUnique, yUnique, regr, corr, theta0, lob, upb);

    testPoints = [X1_vis(:), X2_vis(:)];
    
    predictions = predictor(testPoints, dmodel);
    krigingVis = reshape(predictions, size(X1_vis));
    
    uniquePredictions = predictor(xUnique, dmodel);
    
    % Plot the design and the kriging model predictions.
    plotDesignPoints(X1_vis, X2_vis, krigingVis, xUnique, uniquePredictions, ...
        'Method: Kriging Model');
end
