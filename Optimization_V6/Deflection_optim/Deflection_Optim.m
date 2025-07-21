function [mass, isValid, penalty, W_opt, T_opt, a_opt, b_opt, deflection_mm] = Deflection_Optim(load, H, platform_diam, base_diam, debug, W_init, T_init)
% DEFLECTION_OPTIM Optimizes the SPM linkage design by minimizing mass subject to 
% deflection and stiffness constraints under a specified torque load.
%
%   [mass, isValid, W_opt, T_opt, a_opt, b_opt] = Deflection_Optim(load, H, theta_deg, platform_diam, debug, W_init, T_init)
%
%   Inputs:
%       load         - Applied load (torque in N·m) from Simscape.
%       H            - Mechanism height (m).
%     platform_diam  - Platform diameter (m).
%      base_diam     - Base diameter (m).
%       debug        - (Optional) Boolean flag; if true, prints debug info. Default: false
%       W_init       - (Optional) Initial guess for link width (m). If not provided, a default value is used.
%       T_init       - (Optional) Initial guess for link thickness (m). If not provided, a default value is used.
%
%   Outputs:
%       mass         - Optimized mass (in grams) computed as:
%                      mass = density * (W * T * L) * 1000,
%                      where L is the quarter ellipse arc length.
%       isValid      - Boolean flag: true if the computed deflection is within limits and 
%                      stiffness is within the specified bounds; false otherwise.
%      penalty       - Deflection penalty applied to the mass as a scalar
%       W_opt        - Optimized link width (m)
%       T_opt        - Optimized link thickness (m)
%       a_opt        - Major half-axis of the ellipse (m)
%       b_opt        - Minor half-axis of the ellipse (m)
%     deflection     - Deflection of the linkage (mm)


% Set default value for debug if not provided
if nargin < 5 || isempty(debug)
    debug = false; % Default to no debug output
end


% Material properties (6061 aluminum)
E = 69e9;          % Young's modulus in Pascals
density = 2700;    % Density in kg/m^3

% Calculate base radius based on H and theta (assuming platform_diam represents minor axis b)
%base_rad = H * sind(theta_deg); probably wrong?
base_rad = base_diam/2;
P = load;

% --- Set Bounds for Design Variables (W and T) ---
min_val = 0.001;  % Lower bound in meters
max_val = 0.05;   % Upper bound in meters

% --- Initialize Design Guesses for W and T ---
if nargin < 6 || isempty(W_init)
    %W = min_val + (max_val - min_val) * rand();
    W0 = .01;
else
    W0 = W_init;
end

if nargin < 7 || isempty(T_init)
    %T = min_val + (max_val - min_val) * rand();
    T0 = .005;
else
    T0 = T_init;
end




% %% Test using Optimization Parameters (DOESN'T WORK BECAUSE OPTIMIZATION EXPRESSIONS NOT SUPPORTED IN integralCalc)
% % Create optimization variables for each problem variable
% W = optimvar('W','LowerBound',.002,'UpperBound',.05);
% T = optimvar('T','LowerBound',.001,'UpperBound',.05);
% 
% % Compute the beam deflection
% [deflection_mm, m, ~, ~, L] = curved_beam_deflection(H, platform_diam, base_rad, W, T, load);
% 
% % Create the objective problem
% linprob = optimproblem('Objective', m);
% 
% % Define the deflection constraint
% linprob.Constraints.cons1 = deflection_mm/1000 <= L * .01;
% 
% % Initial Guess
% clear x0
% x0.W = W0; x0.T = T0;
% 
% linsol = solve(linprob, x0);
% 
% % Compute the optimum beam deflection
% W_opt = linsol.W;
% T_opt = linsol.T;
% [deflection_mm, mass, a, b, L] = curved_beam_deflection(H, platform_diam, base_rad, W_opt, T_opt, load);
% deflection_limit = L * .01;
% isValid = deflection_mm/1000 <= deflection_limit;
% if isValid
%     penalty = 1;
% else
%     penalty = 1 + 1e3 * (delfection_mm/1000 - deflection_limit)/deflection_limit;
% end
% 
% W_opt = linsol.W;
% T_opt = linsol.T;
% a_opt = a;
% b_opt = b;

%% Test using fmincon
% Initial guess and bounds for [w, t]
x0 = [W0, T0];
lb = [.003, .0015];
ub = [.05, .025];
options = optimoptions('fmincon','Display','none','Algorithm','sqp');
[x_opt, ~] = fmincon(@obj, x0, [], [], [], [], lb, ub, @nonlcon, options);

% Compute the optimum beam deflection
W_opt = x_opt(1);
T_opt = x_opt(2);
[deflection_mm, m_kg, a, b, Lf] = curved_beam_deflection(H, platform_diam, base_rad, W_opt, T_opt, load);
% deflection_limit = Lf * .01;
deflection_limit = 0.5;  %mm
% deflection_m = deflection_mm .* 1000;
isValid = deflection_mm <= deflection_limit+.01;
if isValid
    penalty = 1.;
else
    penalty = 1 + 1e3 * (deflection_mm - deflection_limit)/deflection_limit;
end
mass = m_kg * 1000;
a_opt = a;
b_opt = b;


% --- Debug Prints Only (no plots) ---
if debug
    fprintf('\n======= OPTIMIZATION RESULTS =======\n');
    fprintf('Optimized Mass: %.2f grams\n', mass);
    fprintf('W (width): %.5f m, T (thickness): %.5f m\n', W_opt, T_opt);
    fprintf('a (major half-axis): %.5f m, b (minor half-axis): %.5f m\n', a_opt, b_opt);
    fprintf('Constraints met: %s\n', mat2str(isValid));
    if ~isValid
        fprintf('Warning: Not all constraints were satisfied\n');
    end
    fprintf('===================================\n\n');
    fprintf('Detailed results:\n');
    fprintf('   W (width)     = %.5f m\n', W_opt);
    fprintf('   T (thickness) = %.5f m\n', T_opt);
    fprintf('   a (major half-axis) = %.5f m\n', a_opt);
    fprintf('   b (minor half-axis) = %.5f m\n', b_opt);
    fprintf('   L (quarter ellipse length) = %.5f m\n', Lf);
    fprintf('   F (force input) = %.2f N\n', load);
    fprintf('   Deflection: %.5f mm (Limit: %.2f mm)\n', deflection_mm, deflection_limit);
    fprintf('   mass = %.3f g\n', mass);
    fprintf('   Penalty = %.3f m\n', penalty);
    % fprintf('Stiffness: %.2f N/m (Required: [%.2f, %.2f] N/m)\n', K, K_min, K_max);
end

% ---------- Objective function (nested) ----------
function m = obj(x)
    w = x(1); t = x(2);
    [~, m, ~, ~, ~] = curved_beam_deflection(H, platform_diam, base_rad, w, t, P);
end

% ---------- Constraint function (nested) ----------
function [c, ceq] = nonlcon(x)
    w = x(1); t = x(2);
    [deflection, ~, ~, ~, ~] = curved_beam_deflection(H, platform_diam, base_rad, w, t, P);
    c1 = deflection - 0.5; % MAKE SURE THIS MATCHES THE DEFLECTION LIMIT ABOVE
    c2 = t/(4*w) - 1; % keep the aspect ratio within reasonable bounds
    c = [c1; c2];
    ceq = [];
end
end

function [deflection_mm, mass, a, b, quarter_circumference] = curved_beam_deflection(h, platform_diam, base_rad, w, t, P)
% curved_beam_deflection calculates the deflection of a quarter-elliptical beam
% using Castigliano's theorem.
%
% Inputs:
% h := height of the SPM (scalar)
% platform_diam := diameter of the platform (scalar) - represents minor axis b
% base_rad := base radius calculated from h and theta (scalar)
% w - cross-section width (m)
% t - cross-section thickness (m)
% P - load applied at the free end (N)
%
% Outputs:
% deflection_castigliano - deflection at the free end (mm)
% mass - mass of the beam (kg)
% a - major half-axis of the ellipse (m)
% b - minor half-axis of the ellipse (m)
% quarter_circumference - 1/4 the circumference of the ellipse (m)

    % Calculate ellipse parameters
    b = platform_diam/2;
    a = sqrt(h^2 + base_rad^2);
    
    % Material properties for 6061 aluminum
    E = 69e9; % Young's modulus (Pa)
    G = E / (2 * (1 + 0.33)); % Shear modulus (Pa)
    rho = 2700; % density of Al in kg/m^3
    
    % Geometric properties
    I = (w * t^3) / 12; % second moment of area for rectangular cross-section
    Kt = w * t^3 / 3; % torsional constant for rectangular cross-section (approximate)
    
    % Calculate 1/4 the circumference of the ellipse (approximation using Ramanujan's formula)
    quarter_circumference = (pi/4) * (3*(a + b) - sqrt((3*a + b)*(a + 3*b)));
    
    % Calculate mass in kg
    mass = rho * (w * t) * quarter_circumference;
    
    % Define numeric functions for y, dy/dx, and ds
    y = @(x) b * sqrt(1 - (x.^2 / a^2));
    dy_dx = @(x) - (b * x) / (a^2 * sqrt(1 - (x.^2 / a^2)));
    ds = @(x) sqrt(1 + dy_dx(x).^2);
    
    % Integrands for bending and torsion components
    Mb = @(x) -P * y(x) .* cos(atan(dy_dx(x))) - P * (a - x) .* sin(atan(dy_dx(x)));
    Mt = @(x) P * y(x) .* sin(atan(dy_dx(x))) - P * (a - x) .* cos(atan(dy_dx(x)));
    
    % Numerical integration
    I1 = integral(@(x) (Mb(x).^2) ./ (E * I) .* ds(x), 0, a);
    I2 = integral(@(x) (Mt(x).^2) ./ (G * Kt) .* ds(x), 0, a);
    
    % Calculate deflection using Castigliano's theorem
    deflection_m = (P / E) * I1 + (P / G) * I2;
    
    % Convert to mm
    deflection_mm_castigliano = deflection_m * 1000;

    % Compute the simple cantilevered beam deflection
    I_cb = w*t^3/12;
    deflection_m_cb = P * quarter_circumference^3 / (3 * E * I_cb);

    % Get the deflection
    deflection_mm = max(deflection_mm_castigliano, deflection_m_cb*1000);


end