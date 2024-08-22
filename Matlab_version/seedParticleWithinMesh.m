function newPos = seedParticleWithinMesh(elements, nodeCoordinates)
    % Randomly select an element from the mesh
    elementIndex = randi(size(elements, 1));
    nodeIndices = elements(elementIndex, :);

    % Get the coordinates of the triangle vertices
    vertices = nodeCoordinates(nodeIndices, :);

    % Generate random barycentric coordinates
    r1 = rand();
    r2 = rand();
    if r1 + r2 > 1
        r1 = 1 - r1;
        r2 = 1 - r2;
    end
    r3 = 1 - r1 - r2;

    % Calculate the new position using the barycentric coordinates
    newPos = r1 * vertices(1, :) + r2 * vertices(2, :) + r3 * vertices(3, :);
end