function [newPos, newVel] = applyReflectiveBoundary(pos, vel, minLat, maxLat, minLon, maxLon)
    newPos = pos;
    newVel = vel;

    % Reflect latitude
    if newPos(2) < minLat
        newPos(2) = 2*minLat - newPos(2);
        newVel(2) = -newVel(2);
    elseif newPos(2) > maxLat
        newPos(2) = 2*maxLat - newPos(2);
        newVel(2) = -newVel(2);
    end

    % Reflect longitude
    if newPos(1) < minLon
        newPos(1) = 2*minLon - newPos(1);
        newVel(1) = -newVel(1);
    elseif newPos(1) > maxLon
        newPos(1) = 2*maxLon - newPos(1);
        newVel(1) = -newVel(1);
    end
end