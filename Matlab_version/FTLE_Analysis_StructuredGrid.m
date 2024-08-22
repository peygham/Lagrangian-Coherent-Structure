%% LCS 
% This MATLAB code is designed to compute the Finite-Time Lyapunov Exponent (FTLE) fields, 
% both forward and backward, over a structured grid. FTLE fields are used to identify 
% Lagrangian Coherent Structures (LCS) in fluid flows, which can represent regions of 
% high stretching or compression. These structures are often interpreted as repelling 
% or attracting manifolds that play a crucial role in the transport and mixing of particles in the flow.

% The boundary condition implemented in this code is a reflective boundary condition
% Integration Approach: Runge-Kutta 4th Order (RK4)
% pgh@akvaplan.niva.no

clear; clc; close all;
filePath = 'C:\RFF_NordLand\LCS_PGH\LCS4fv\data_March';
load MarMesh

% S.X, S.Y, and S.time are correctly defined in MarMesh
lon_grid = S.X;
lat_grid = S.Y;
time = S.time;

% Pre-allocate grid data
num_days = 1; % Number of days to process
u_grid = [];
v_grid = [];

% Load velocity fields and concatenate them
for day_num = 23:23 % Example: day 23
    disp(['Day number = ', num2str(day_num)]);
    dataFile = fullfile(filePath, ['Velocity' num2str(day_num, '%02d') '.mat']);
    load(dataFile, 'field');
    u_grid = cat(1, u_grid, field.u);
    v_grid = cat(1, v_grid, field.v);
    clear field;
end

% Create a mask where velocity data is valid
valid_mask = squeeze(~isnan(u_grid(1,:,:)) & ~isnan(v_grid(1,:,:)));

% Initial particle positions
x0 = lon_grid;
y0 = lat_grid;

%% Parameters for FTLE calculation
dt = 1/(60*60); % Time step for integration (1 second in hour units)
integration_time = 24; % Total integration time (hours)
num_steps = round(integration_time / dt);

%% FTLE FORWARD
% Pre-allocate FTLE forward field
ftle_forward = nan(size(lon_grid));

% Precompute constants for the reflective boundary conditions
lon_min = min(lon_grid(:));
lon_max = max(lon_grid(:));
lat_min = min(lat_grid(:));
lat_max = max(lat_grid(:));
lon_diff = lon_max - lon_min;
lat_diff = lat_max - lat_min;

% Forward integration using RK4 method
for t = 1:num_steps
    if mod(t, 10) == 0
        disp(['Forward integration step: ', num2str(t)]);
    end
    
    % Current velocity fields
    u_curr = squeeze(u_grid(mod(t-1, size(u_grid, 1)) + 1, :, :));
    v_curr = squeeze(v_grid(mod(t-1, size(v_grid, 1)) + 1, :, :));
    
    % Replace NaNs with zero for interpolation safety
    u_curr(isnan(u_curr)) = 0;
    v_curr(isnan(v_curr)) = 0;
    
    % Perform vectorized RK4 step
    [dx, dy] = rk4_step(x0, y0, u_curr, v_curr, dt, lon_grid, lat_grid);
    
    % Update positions
    x_new = x0 + dx;
    y_new = y0 + dy;

    % Reflective boundary conditions handling
    x_new = mod(x_new - lon_min, lon_diff) + lon_min;
    y_new = mod(y_new - lat_min, lat_diff) + lat_min;

    % Ensure indices are within valid range
    x_idx = round((x_new - lon_min) / lon_diff * (size(valid_mask, 2) - 1)) + 1;
    y_idx = round((y_new - lat_min) / lat_diff * (size(valid_mask, 1) - 1)) + 1;
    
    x_idx = min(max(x_idx, 1), size(valid_mask, 2));
    y_idx = min(max(y_idx, 1), size(valid_mask, 1));

    % Create a logical mask for valid positions
    valid_moves = valid_mask(sub2ind(size(valid_mask), y_idx(:), x_idx(:)));

    % Update positions only if valid
    x0(valid_moves) = x_new(valid_moves);
    y0(valid_moves) = y_new(valid_moves);
end

% Compute gradients in a vectorized manner
[grad_x_x0, grad_y_x0] = gradient(x0, lon_grid(1,:), lat_grid(:,1));
[grad_x_y0, grad_y_y0] = gradient(y0, lon_grid(1,:), lat_grid(:,1));

% Compute FTLE forward in a vectorized manner
F11 = grad_x_x0; F12 = grad_y_x0;
F21 = grad_x_y0; F22 = grad_y_y0;
C11 = F11.^2 + F21.^2; C12 = F11.*F12 + F21.*F22;
C22 = F12.^2 + F22.^2;
lambda_max = 0.5 * (C11 + C22 + sqrt((C11 - C22).^2 + 4*C12.^2));
ftle_forward = log(sqrt(lambda_max)) / (integration_time * 3600);

%% FTLE BACKWARD
% Pre-allocate FTLE backward field
ftle_backward = nan(size(lon_grid));

% Reverse the velocity fields for backward integration
u_grid_reverse = -u_grid;
v_grid_reverse = -v_grid;

% Initialize particle positions for backward integration
x0 = lon_grid;
y0 = lat_grid;

% Backward integration using RK4 method
for t = 1:num_steps
    if mod(t, 10) == 0
        disp(['Backward integration step: ', num2str(t)]);
    end
    
    u_curr = squeeze(u_grid_reverse(mod(t-1, size(u_grid_reverse, 1)) + 1, :, :));
    v_curr = squeeze(v_grid_reverse(mod(t-1, size(v_grid_reverse, 1)) + 1, :, :));
    
    u_curr(isnan(u_curr)) = 0;
    v_curr(isnan(v_curr)) = 0;
    
    % Perform vectorized RK4 step
    [dx, dy] = rk4_step(x0, y0, u_curr, v_curr, dt, lon_grid, lat_grid);
    
    x_new = x0 + dx;
    y_new = y0 + dy;

    % Reflective boundary conditions handling
    x_new = mod(x_new - lon_min, lon_diff) + lon_min;
    y_new = mod(y_new - lat_min, lat_diff) + lat_min;

    % Ensure indices are within valid range
    x_idx = round((x_new - lon_min) / lon_diff * (size(valid_mask, 2) - 1)) + 1;
    y_idx = round((y_new - lat_min) / lat_diff * (size(valid_mask, 1) - 1)) + 1;
    
    x_idx = min(max(x_idx, 1), size(valid_mask, 2));
    y_idx = min(max(y_idx, 1), size(valid_mask, 1));

    % Create a logical mask for valid positions
    valid_moves = valid_mask(sub2ind(size(valid_mask), y_idx(:), x_idx(:)));

    % Update positions only if valid
    x0(valid_moves) = x_new(valid_moves);
    y0(valid_moves) = y_new(valid_moves);
end

% Compute gradients in a vectorized manner
[grad_x_x0_b, grad_y_x0_b] = gradient(x0, lon_grid(1,:), lat_grid(:,1));
[grad_x_y0_b, grad_y_y0_b] = gradient(y0, lon_grid(1,:), lat_grid(:,1));

% Compute FTLE backward in a vectorized manner
F11_b = grad_x_x0_b; F12_b = grad_y_x0_b;
F21_b = grad_x_y0_b; F22_b = grad_y_y0_b;
C11_b = F11_b.^2 + F21_b.^2; C12_b = F11_b.*F12_b + F21_b.*F22_b;
C22_b = F12_b.^2 + F22_b.^2;
lambda_max_b = 0.5 * (C11_b + C22_b + sqrt((C11_b - C22_b).^2 + 4*C12_b.^2));
ftle_backward = log(sqrt(lambda_max_b)) / (integration_time * 3600);

%% Plot FTLE
% Plotting Forward and Backward FTLE Fields
ftleF = ftle_forward;
ftleF(~valid_mask) = nan;
figure;
pcolor(lon_grid, lat_grid, ftleF); 
shading interp;
colorbar;
cmocean('delta');
title('Forward FTLE (Repelling Structures)');
xlabel('Longitude');
ylabel('Latitude');

ftleB = ftle_backward;
ftleB(~valid_mask) = nan;
figure;
pcolor(lon_grid, lat_grid, ftleB);
shading interp;
colorbar;
cmocean('delta');
title('Backward FTLE (Attracting Structures)');
xlabel('Longitude');
ylabel('Latitude');