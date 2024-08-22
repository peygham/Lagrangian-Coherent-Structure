# FastFTLE_LCS_Analysis with Parallel Processing
# This MATLAB code is designed to compute the Finite-Time Lyapunov Exponent (FTLE) fields, 
# both forward and backward, over a structured grid. FTLE fields are used to identify 
# Lagrangian Coherent Structures (LCS) in fluid flows, which can represent regions of 
# high stretching or compression. These structures are often interpreted as repelling 
# or attracting manifolds that play a crucial role in the transport and mixing of particles in the flow.

# The boundary condition implemented in this code is NOT a reflective boundary condition
# Integration Approach: Runge-Kutta 4th Order (RK4)
# pgh@akvaplan.niva.no
# Aug.2024

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import map_coordinates, sobel
import h5py

# Load the grid data
with h5py.File('./input/MarMesh23.mat', 'r') as f:
    lon_grid = np.array(f['S']['X'][:])  # Ensure data is loaded as a copy
    lat_grid = np.array(f['S']['Y'][:])  # Ensure data is loaded as a copy
    time = np.array(f['S']['time'][:])   # Ensure data is loaded as a copy

# Convert grid arrays to the correct type (if necessary)
lon_grid = np.array(lon_grid, dtype=np.float32)
lat_grid = np.array(lat_grid, dtype=np.float32)

# Ensure that lat_grid and lon_grid are strictly ascending
if np.any(np.diff(lat_grid[:, 0]) <= 0):
    sorted_indices_lat = np.argsort(lat_grid[:, 0])
    lat_grid = lat_grid[sorted_indices_lat, :]

if np.any(np.diff(lon_grid[0, :]) <= 0):
    sorted_indices_lon = np.argsort(lon_grid[0, :])
    lon_grid = lon_grid[:, sorted_indices_lon]

# Create a mask for valid points (where there are no NaNs)
valid_mask = ~np.isnan(lon_grid) & ~np.isnan(lat_grid)

# Parameters for FTLE calculation
dt = 1/(60)  # Time step for integration (1 min in hour units)
integration_time = 1  # Total integration time for each chunk (1 hour)
num_steps = int(integration_time / dt)

# Precompute constants for the reflective boundary conditions
lon_min, lon_max = lon_grid.min(), lon_grid.max()
lat_min, lat_max = lat_grid.min(), lat_grid.max()
lon_diff = lon_max - lon_min
lat_diff = lat_max - lat_min

# FTLE fields
ftle_forward = np.full(lon_grid.shape, np.nan)
ftle_backward = np.full(lon_grid.shape, np.nan)

def rk4_step(x, y, u, v, dt, lon_grid, lat_grid, valid_mask):
    def interpolate(data, x, y):
        xi = np.array([(y - lat_grid.min()) / (lat_grid.max() - lat_grid.min()) * (data.shape[0] - 1),
                       (x - lon_grid.min()) / (lon_grid.max() - lon_grid.min()) * (data.shape[1] - 1)])
        return map_coordinates(data, xi, order=1, mode='nearest')

    k1_x = interpolate(u, x, y)
    k1_y = interpolate(v, x, y)
    
    k2_x = interpolate(u, x + 0.5 * dt * k1_x, y + 0.5 * dt * k1_y)
    k2_y = interpolate(v, x + 0.5 * dt * k1_x, y + 0.5 * dt * k1_y)
    
    k3_x = interpolate(u, x + 0.5 * dt * k2_x, y + 0.5 * dt * k2_y)
    k3_y = interpolate(v, x + 0.5 * dt * k2_x, y + 0.5 * dt * k2_y)
    
    k4_x = interpolate(u, x + dt * k3_x, y + dt * k3_y)
    k4_y = interpolate(v, x + dt * k3_x, y + dt * k3_y)
    
    dx = (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6.0 * dt
    dy = (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6.0 * dt
    
    # Apply the valid mask to avoid NaN areas
    dx[~valid_mask] = np.nan
    dy[~valid_mask] = np.nan
    
    return dx, dy

# Process velocity fields in hourly chunks
for hour in range(1, 2):  # Loop through 24 hours
    print(f'Processing hour: {hour}')
    
    with h5py.File('./input/VelocityMar5_U_23.mat', 'r') as f_u, h5py.File('./input/VelocityMar5_V_23.mat', 'r') as f_v:
        # Load only the necessary chunk of data for the current hour
        u_chunk = f_u['V'][:, :, hour-1]
        v_chunk = f_v['V'][:, :, hour-1]
    
        # Reverse the velocity fields for backward integration
        u_chunk_reverse = -u_chunk
        v_chunk_reverse = -v_chunk

        # Initialize particle positions for forward and backward integration
        x0 = lon_grid.copy()
        y0 = lat_grid.copy()
        x0_b = lon_grid.copy()
        y0_b = lat_grid.copy()

        for t in range(num_steps):
            if t % 10 == 0:
                print(f'Integration step: {t + (hour-1)}')

            # FORWARD INTEGRATION
            u_curr = u_chunk
            v_curr = v_chunk

            dx, dy = rk4_step(x0, y0, u_curr, v_curr, dt, lon_grid, lat_grid, valid_mask)

            x_new = np.mod(x0 + dx - lon_min, lon_diff) + lon_min
            y_new = np.mod(y0 + dy - lat_min, lat_diff) + lat_min

            x0, y0 = x_new, y_new

            # BACKWARD INTEGRATION
            u_curr_reverse = u_chunk_reverse
            v_curr_reverse = v_chunk_reverse

            dx_b, dy_b = rk4_step(x0_b, y0_b, u_curr_reverse, v_curr_reverse, dt, lon_grid, lat_grid, valid_mask)

            x_new_b = np.mod(x0_b + dx_b - lon_min, lon_diff) + lon_min
            y_new_b = np.mod(y0_b + dy_b - lat_min, lat_diff) + lat_min

            x0_b, y0_b = x_new_b, y_new_b

    # Compute gradients for FTLE after processing each chunk
    grad_x_x0 = sobel(x0, axis=1)
    grad_y_x0 = sobel(x0, axis=0)
    grad_x_y0 = sobel(y0, axis=1)
    grad_y_y0 = sobel(y0, axis=0)

    F11 = grad_x_x0
    F12 = grad_y_x0
    F21 = grad_x_y0
    F22 = grad_y_y0

    C11 = F11**2 + F21**2
    C12 = F11 * F12 + F21 * F22
    C22 = F12**2 + F22**2

    lambda_max = 0.5 * (C11 + C22 + np.sqrt((C11 - C22)**2 + 4 * C12**2))
    ftle_forward = np.log(np.sqrt(lambda_max)) / (integration_time * 3600)

    # Compute gradients for backward FTLE
    grad_x_x0_b = sobel(x0_b, axis=1)
    grad_y_x0_b = sobel(x0_b, axis=0)
    grad_x_y0_b = sobel(y0_b, axis=1)
    grad_y_y0_b = sobel(y0_b, axis=0)

    F11_b = grad_x_x0_b
    F12_b = grad_y_x0_b
    F21_b = grad_x_y0_b
    F22_b = grad_y_y0_b

    C11_b = F11_b**2 + F21_b**2
    C12_b = F11_b * F12_b + F21_b * F22_b
    C22_b = F12_b**2 + F22_b**2

    lambda_max_b = 0.5 * (C11_b + C22_b + np.sqrt((C11_b - C22_b)**2 + 4 * C12_b**2))
    ftle_backward = np.log(np.sqrt(lambda_max_b)) / (integration_time * 3600)

    # Plot FTLE for each chunk
    plt.figure()
    plt.pcolor(lon_grid, lat_grid, ftle_forward)
    plt.colorbar()
    plt.title(f'Forward FTLE (Repelling Structures) - Hours {hour}')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.savefig(f'Forward_FTLE_Hour_{hour}.eps', format='eps')

    plt.figure()
    plt.pcolor(lon_grid, lat_grid, ftle_backward)
    plt.colorbar()
    plt.title(f'Backward FTLE (Attracting Structures) - Hours {hour}')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.savefig(f'Backward_FTLE_Hour_{hour}.eps', format='eps')

    # Reset particle positions for next chunk
    x0 = lon_grid.copy()
    y0 = lat_grid.copy()
    x0_b = lon_grid.copy()
    y0_b = lat_grid.copy()

