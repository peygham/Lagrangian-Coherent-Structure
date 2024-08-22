function baryCoords = calculateBarycentricCoordinates(pt, tri)
    % Implement barycentric coordinate calculation
    x = pt(1); y = pt(2);
    x1 = tri(1,1); y1 = tri(1,2);
    x2 = tri(2,1); y2 = tri(2,2);
    x3 = tri(3,1); y3 = tri(3,2);
    
    denominator = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3);
    a = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / denominator;
    b = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / denominator;
    baryCoords = [a, b, 1 - a - b];
end