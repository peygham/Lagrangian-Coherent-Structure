function valid = isValidPosition(position, nodeCoordinates, elements)
    % Initialize validity to false
    valid = false;
    
    % Loop through all elements
    for i = 1:size(elements, 1)
        % Get the nodes that make up this element
        nodeIndices = elements(i, :);
        
        % Get the coordinates of the nodes
        vertices = nodeCoordinates(nodeIndices, :);
        
        % Check if the point is inside the triangle using the barycentric method
        A = [vertices(1, :) - position; vertices(2, :) - position; vertices(3, :) - position];
        T = A(1:2, :) \ (vertices(3, :) - position)';
        
        % The point is inside the triangle if all the barycentric coordinates are between 0 and 1
        if all(T >= 0) && all(T <= 1) && sum(T) <= 1
            valid = true;
            return;
        end
    end
end