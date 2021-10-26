function perimeter = normalized_L(XY)
    % Compute the perimeter of an input curve represented by an N x 2 matrix
    segments = circshift(XY, -1) - XY; % Segment of the polygon
    perimeter = sum(sqrt(sum(segments .* segments, 2)));  % Sum of segment lengths
end
