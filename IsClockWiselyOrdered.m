function isClockWise = IsClockWiselyOrdered(points) % https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order/18472899
    lenPoints = length(points);
    isClockWise = false;
    signedArea = 0.0; % See https://stackoverflow.com/questions/451426/how-do-i-calculate-the-area-of-a-2d-polygon
    for i=1:2:lenPoints
        A = [points(i) points(i+1)];
        if i==(lenPoints-1)
            B = [points(1) points(2)];
        else
            B = [points(i+2) points(i+3)];
        end
        signedArea = signedArea + det([A; B]);% See http://mathworld.wolfram.com/PolygonArea.html
    end % signedArea = signedArea*0.5; % Notice the area must be divided by 2 to be correct! But as we do not need it...
    if signedArea<0 % If it is negative the points are in clockwise order
        isClockWise = true;
    end
end
