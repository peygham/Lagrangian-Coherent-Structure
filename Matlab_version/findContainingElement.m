function [elementIndex, baryCoords] = findContainingElement(position, nodeCoordinates, elements)
    numElements = size(elements, 1);
    elementIndex = NaN;  % Default value if no containing element is found
    baryCoords = NaN(1, 3);  % Default barycentric coordinates

    for i = 1:numElements
        % Get the nodes of the current element
        nodeIndices = elements(i, :);
        nodes = nodeCoordinates(nodeIndices, :);

        % Check if the point is inside the triangle
        if isPointInTriangle(position, nodes)
            % Calculate barycentric coordinates
            baryCoords = calculateBarycentricCoordinates(position, nodes);
            elementIndex = i;
            return;
        end
    end
end
