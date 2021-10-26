function [nsdvec] = GetShapeVector(x_points, y_points, matchOrder)
%GETSHAPEVECTOR Function to compute the Normalized Shape-Descriptor Vector
%   This function computes the normalized shape-descriptor vector (NDSV) of
%   a given trajectory. The trajectory must be described by 2 column-vecs.
%   The first one has to contains the x's coordinates [x1; x2...] while the
%   second one must contain the y's coordinates [y1; y2...]. They must be 
%   ordered (no matter the sense)  The output is a vector with [1, angle 
%   with the next segment / 360º...,length of the next segment divided by 
%   that of the longest... If the path is closed, OK. Otrw., the final 
%   segment is virtually considered
%   Do matchOrder = 1 to ensure that the final vector has a comparable order
    points = reshape([x_points, y_points]', 1, 2*length(x_points)); % NUEVO
    %points Si, da como si pusieras uno tras otro
    lenPoints = length(points);
    if size(points, 1)~= 1 || mod(lenPoints, 2) ~= 0
        error('The input must be a list of pairs of coordinates');
    end

    if IsClockWiselyOrdered(points)
        points = ReversePolygon(points, lenPoints);
    end
    [longest_line, L] = SeekLongestSide(points, lenPoints);
    
    nsdvec = ComputeNSDVector(points, lenPoints, longest_line, L);
    if exist('matchOrder', 'var') && matchOrder == true
        nsdvec = MatchOrder(nsdvec, lenPoints);
    end
end

% Internal auxiliary functions

function points = ReversePolygon(points, lenPoints) % Iterating in inverse order is more efficient but more error-prone too
    focus = 1;
    for i=1:1:lenPoints/4
        swapX = points(focus);
        swapY = points(focus+1);
        points(focus) = points(lenPoints-focus);
        points(focus+1) = points(lenPoints-focus+1);
        points(lenPoints-focus) = swapX;
        points(lenPoints-focus+1) = swapY;
        focus = focus + 2;
    end
end

function [index, max_len] = SeekLongestSide(points, lenPoints)
    index = 1;
    max_len = -1;
    for i=1:2:lenPoints
        A = [points(i), points(i+1)];
        if i<(lenPoints-1)
            B = [points(i+2), points(i+3)];
        else
            B = [points(1), points(2)];
        end
        dist = norm(A-B);
        if dist > max_len
            max_len = dist;
            index = i;
        end
    end
end

function [nsdvec] = ComputeNSDVector(points, lenPoints, longest_line, L)
    nsdvec = zeros(1, lenPoints);
    focus = longest_line;
    
    A = [points(focus) points(focus+1)];
    if focus<(lenPoints-1)
        B = [points(focus+2), points(focus+3)];
    else
        B = [points(1), points(2)];
    end
    
    for i=1:1:lenPoints
        if mod(i, 2)==1
            F = A; % For concavity detection of every single angle
            
            AB = (B-A)/norm(B-A);
            nsdvec(i) = norm(B-A)/L;
            focus = focus + 2;
            if focus > (lenPoints - 1)
                focus = 1;
            end
            A = B; % [points(focus) points(focus+1)]; % Moving...
            if focus<(lenPoints-1)
                B = [points(focus+2), points(focus+3)];
            else
                B = [points(1), points(2)];
            end
        else
            CD = (A-B)/norm(A-B); % Force the same origin
            angle = acosd(AB*(CD'));
            if CheckLocalConcavity(F, A, B) % Concavity check (A should be the former G and B should be the former H as a result of moving)
                angle = 360 - angle;
            end
            nsdvec(i) = angle/360;
        end
    end
end

function isConcave=CheckLocalConcavity(F, G, H) % See: https://en.wikipedia.org/wiki/Curve_orientation
    isConcave = false;                          % See: https://matlabgeeks.com/tips-tutorials/computational-geometry/check-convexity-of-polygon/
    O = [1 F; 1 G; 1 H];                        % See: https://community.khronos.org/t/algorithm-to-detect-the-type-of-a-polygon/55974/4
    detO = det(O); % As we force positively oriented polygons (counterclockwise)
    if detO<0 % If determinant of orientation matrix for local points is negative
       isConcave = true; % Concave sequence of points
    end                                         % See: https://en.wikipedia.org/wiki/Determinant#Orientation_of_a_basis
end

function nsdvec = MatchOrder(nsdvec, lenNSDVec) % Order the NSDVector so that it is always the same for polygons even though they have more than a longest side
    potentialInits = []; % Of course, this info could be recorded while computing the NSD vector itself, but it is done here for modularity
    for i=1:2:length(nsdvec)
        if(nsdvec(i)==1.0)
           potentialInits = [potentialInits i];
        end
    end
    numInits = length(potentialInits);
    if numInits>1 % Act only if necessary, otherwise, the vector starts with 1.0 for its longest side
        focuses = AdvanceFocuses(potentialInits, lenNSDVec); % We start after moving the start, as they all should be 1 
        steps = 1;
        while numInits>1 && steps<lenNSDVec
            participants = nsdvec(focuses);
            current_reference = max(participants);
            id_alive = find(participants==current_reference);
            numInits = length(id_alive);
            potentialInits = potentialInits(id_alive);
            focuses = AdvanceFocuses(focuses(id_alive), lenNSDVec);
            steps = steps+1;
        end
        if ~isempty(potentialInits) % Hotfix to avoid that NaN's linked to degenerated shapes leave no point to start at.
            nsdvec = RefactorVector(nsdvec, potentialInits(1), lenNSDVec);
        end
    end
end

function focuses = AdvanceFocuses(focuses, lenNSDVec)
    for i=1:length(focuses)
        focuses(i) = focuses(i)+1;
        if focuses(i)>lenNSDVec
            focuses(i) = 1;
        end
    end
end

function finalVec = RefactorVector(nsdvec, focus, lenNSDVec)
    finalVec = zeros(1, lenNSDVec);
    for i=1:1:length(nsdvec)
        finalVec(i) = nsdvec(focus);
        focus = focus + 1;
        if focus>lenNSDVec
            focus = 1;
        end
    end
end
