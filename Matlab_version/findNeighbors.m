function neighbors = findNeighbors(elements, index)
    % Find the coordinates of all nodes in the elements
    x = elements(:, 1);
    y = elements(:, 2);

    % Check if the index is within the valid range
    if index < 1 || index > length(x)
        error('Index exceeds the bounds of the coordinate arrays.');
    end

    % Calculate the distances from the point at 'index' to all other points
    dist = sqrt((x - x(index)).^2 + (y - y(index)).^2);

    % Define a distance threshold for neighbors (e.g., based on grid resolution)
    distanceThreshold = 1e-3;  % You may need to adjust this threshold

    % Find all points that are within the distance threshold
    neighbors = find(dist < distanceThreshold & dist > 0); % Exclude the point itself
end