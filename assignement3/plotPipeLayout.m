function plotPipeLayout(L_e, D_p, N)
% plotSerpentinePipe plots a top–down view of a serpentine pipe.
%
%   plotSerpentinePipe(L_e, D_p, N) plots a pipe of total length L_e,
%   pipe “diameter” D_p (which also defines the offset between the two edges),
%   and N 180° bends. For this example the pipe is constructed as a series 
%   of horizontal segments (with equal lengths) connected by half–circle arcs.
%
%   The arcs alternate in orientation: odd–numbered bends curve upward
%   (using radius R_outer = 1.5*D_p) and even–numbered bends curve downward
%   (using radius R_inner = 0.5*D_p). The “top line” (upper edge of the pipe)
%   is obtained by offsetting the computed (bottom) path by D_p in the y–direction.
%
%   NOTE: This implementation is one interpretation of the design. The
%   allocation of straight–and–arc lengths is chosen so that the sum of the
%   horizontal (straight) segments is L_e minus the sum of the arc lengths.
%

% Define bend radii
R_outer = 1.5 * D_p; % for odd bends (upward)
R_inner = 0.5 * D_p; % for even bends (downward)

% Compute arc lengths for each bend (a half circle = pi*R)
arcLengths = zeros(N,1);
for i = 1:N
    if mod(i,2)==1
        arcLengths(i) = pi * R_outer;
    else
        arcLengths(i) = pi * R_inner;
    end
end
totalArcLength = sum(arcLengths);

% The remaining length is allocated to horizontal (straight) segments.
% There are N+1 straight segments.
L_horiz = L_e - totalArcLength;
if L_horiz < 0
    error('L_e is too short for the chosen number of bends and D_p.');
end
L_seg = L_horiz / (N+1);

% Initialize arrays to store the coordinates of the bottom (center) path
x_bottom = [];
y_bottom = [];

% Start at the origin, going east (right)
current_point = [0, 0];
current_dir = [1, 0];  % unit vector pointing east

% Add initial point and first straight segment
x_bottom(end+1) = current_point(1);
y_bottom(end+1) = current_point(2);
current_point = current_point + L_seg * current_dir;
x_bottom(end+1) = current_point(1);
y_bottom(end+1) = current_point(2);

% Loop over each bend to build the serpentine path
for i = 1:N
    if mod(i,2)==1
        % --- Odd bend: upward U-turn ---
        % For a U-turn that changes direction from east to west with horizontal
        % tangents, we choose the following construction:
        %
        % Let the arc start at A = current_point and end at B.
        % For an upward bend we want B = A + [0, 2*R_outer] so that the chord is vertical.
        % The circle that contains A and B with horizontal tangents at both ends
        % has center C = A + [R_outer, R_outer].
        %
        % Parameterize theta from -3*pi/4 to 3*pi/4.
        R = R_outer;
        C = current_point + [R, R];
        theta = linspace(-3*pi/4, 3*pi/4, 50);
        arc_x = C(1) + R*cos(theta);
        arc_y = C(2) + R*sin(theta);
        % Update current point to end of arc
        current_point = [arc_x(end), arc_y(end)];
        % After an upward U-turn the new direction is to the left (west)
        current_dir = [-1, 0];
        
    else
        % --- Even bend: downward U-turn ---
        % For a downward bend we mirror the above.
        % Let the arc start at A = current_point and end at B = A + [0, -2*R_inner].
        % The circle will have center C = A + [R_inner, -R_inner].
        % Parameterize theta from 3*pi/4 down to -3*pi/4.
        R = R_inner;
        C = current_point + [R, -R];
        theta = linspace(3*pi/4, -3*pi/4, 50);
        arc_x = C(1) + R*cos(theta);
        arc_y = C(2) + R*sin(theta);
        % Update current point to end of arc
        current_point = [arc_x(end), arc_y(end)];
        % After a downward U-turn the direction remains to the left (west)
        current_dir = [-1, 0];
    end
    
    % Append the arc points to the bottom path
    x_bottom = [x_bottom, arc_x];
    y_bottom = [y_bottom, arc_y];
    
    % For all but the last bend, add the next straight segment.
    if i < N
        current_point = current_point + L_seg * current_dir;
        x_bottom(end+1) = current_point(1);
        y_bottom(end+1) = current_point(2);
    end
end

% Add the final horizontal segment after the last bend
current_point = current_point + L_seg * current_dir;
x_bottom(end+1) = current_point(1);
y_bottom(end+1) = current_point(2);

% Now compute the top line as an offset of the bottom line.
% (Here we simply add D_p in the y–direction to every point.)
x_top = x_bottom;
y_top = y_bottom + D_p;

% Plot the resulting pipe edges.
figure; hold on; axis equal;
plot(x_bottom, y_bottom, 'b-', 'LineWidth', 2);
plot(x_top, y_top, 'r-', 'LineWidth', 2);
title('Serpentine Pipe');
xlabel('X'); ylabel('Y');
legend('Bottom Edge (Centerline)','Top Edge (Offset)');
grid on;

end
