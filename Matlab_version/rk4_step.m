function [dx, dy] = rk4_step(x, y, u_curr, v_curr, dt, lon_grid, lat_grid)
    % Ensure the interpolated velocity fields are valid
    k1_x = interp2(lon_grid, lat_grid, u_curr, x, y, 'linear', NaN);
    k1_y = interp2(lon_grid, lat_grid, v_curr, x, y, 'linear', NaN);
    
    k2_x = interp2(lon_grid, lat_grid, u_curr, x + 0.5*dt*k1_x, y + 0.5*dt*k1_y, 'linear', NaN);
    k2_y = interp2(lon_grid, lat_grid, v_curr, x + 0.5*dt*k1_x, y + 0.5*dt*k1_y, 'linear', NaN);
    
    k3_x = interp2(lon_grid, lat_grid, u_curr, x + 0.5*dt*k2_x, y + 0.5*dt*k2_y, 'linear', NaN);
    k3_y = interp2(lon_grid, lat_grid, v_curr, x + 0.5*dt*k2_x, y + 0.5*dt*k2_y, 'linear', NaN);
    
    k4_x = interp2(lon_grid, lat_grid, u_curr, x + dt*k3_x, y + dt*k3_y, 'linear', NaN);
    k4_y = interp2(lon_grid, lat_grid, v_curr, x + dt*k3_x, y + dt*k3_y, 'linear', NaN);
    
    dx = dt / 6 * (k1_x + 2*k2_x + 2*k3_x + k4_x);
    dy = dt / 6 * (k1_y + 2*k2_y + 2*k3_y + k4_y);
    
    % Handle NaN values by setting dx and dy to NaN where interpolation failed
    dx(isnan(k1_x) | isnan(k2_x) | isnan(k3_x) | isnan(k4_x)) = NaN;
    dy(isnan(k1_y) | isnan(k2_y) | isnan(k3_y) | isnan(k4_y)) = NaN;
end