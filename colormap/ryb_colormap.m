function cmap = ryb_colormap(n)
%RYB_COLORMAP Creates a red-yellow-blue colormap with variable step size
%
%   CMAP = RYB_COLORMAP(N) generates a colormap with N colors transitioning
%   smoothly from red through yellow to blue. The colors are designed to be
%   nicely distinguishable with good perceptual spacing.
%
%   Input:
%       n - Number of colors in the colormap (step size). Must be >= 2.
%           Default is 64 if not specified.
%
%   Output:
%       cmap - Nx3 array of RGB values (range 0-1)
%
%   Example:
%       colormap(ryb_colormap(64))
%       imagesc(rand(10,10))
%       colorbar
%
%       % Use with different step sizes
%       colormap(ryb_colormap(128))  % Higher resolution
%       colormap(ryb_colormap(32))   % Lower resolution

if nargin < 1
    n = 64; % Default to 64 colors
end

if n < 2
    error('Number of colors must be at least 2');
end

% Create parameter vector from 0 to 1
t = linspace(0, 1, n)';

% Use perceptually better spacing - emphasize the yellow region slightly
% This helps with distinguishability in the middle range
t_perceptual = t.^0.85;

% Define key color waypoints for smooth, distinguishable transitions
% Red -> Orange -> Yellow -> Green -> Cyan -> Blue
waypoints = [
    1.0, 0.0, 0.0;  % Red
    1.0, 0.5, 0.0;  % Orange
    1.0, 1.0, 0.0;  % Yellow
    0.5, 1.0, 0.3;  % Yellow-green
    0.0, 0.8, 0.6;  % Cyan-green
    0.0, 0.4, 1.0;  % Blue-cyan
    0.0, 0.0, 1.0   % Blue
];

% Normalize waypoint positions (0 to 1)
waypoint_pos = linspace(0, 1, size(waypoints, 1))';

% Initialize colormap
cmap = zeros(n, 3);

% Interpolate between waypoints
for i = 1:n
    % Find which segment we're in
    pos = t_perceptual(i);
    
    % Find the two waypoints to interpolate between
    idx = find(waypoint_pos <= pos, 1, 'last');
    
    if isempty(idx) || idx == length(waypoint_pos)
        if idx == length(waypoint_pos)
            cmap(i, :) = waypoints(end, :);
        else
            cmap(i, :) = waypoints(1, :);
        end
    else
        % Interpolate between waypoints idx and idx+1
        t_local = (pos - waypoint_pos(idx)) / (waypoint_pos(idx+1) - waypoint_pos(idx));
        % Use smooth interpolation (smoothstep-like)
        t_smooth = t_local^2 * (3 - 2*t_local);
        cmap(i, :) = (1 - t_smooth) * waypoints(idx, :) + t_smooth * waypoints(idx+1, :);
    end
end

% Ensure values are in valid range [0, 1]
cmap = max(0, min(1, cmap));

end

