function [x_valid, y_valid] = nearest_valid(x, y, valid_lon, valid_lat)
    % Initialize output arrays
    x_valid = zeros(size(x));
    y_valid = zeros(size(y));
    
    % Iterate over each out-of-bounds point
    for i = 1:numel(x)
        % Compute the distances to all valid points
        dist_squared = (valid_lon - x(i)).^2 + (valid_lat - y(i)).^2;
        
        % Find the index of the minimum distance
        [~, closest_idx] = min(dist_squared);
        
        % Set the closest valid point
        x_valid(i) = valid_lon(closest_idx);
        y_valid(i) = valid_lat(closest_idx);
    end
end