% Function to perform farthest point sampling
function idx = farthestPointSampling(ptCloud, numPoints)
    % Extract point locations
    points = ptCloud.Location;
    numOriginalPoints = size(points, 1);
    
    if numPoints > numOriginalPoints
        error('Number of points to sample exceeds the total number of points in the point cloud.');
    end
    
    % Initialize variables
    idx = zeros(numPoints, 1);
    selectedIdx = 1; % Start with the first point
    idx(1) = selectedIdx;
    
    % Create a matrix to store distances to the nearest selected point
    distMatrix = inf(numOriginalPoints, 1);
    
    % Iterate to find farthest points
    for i = 2:numPoints
        % Compute distances from all points to the nearest selected point
        distToNearestSelected = sqrt(sum((points - points(selectedIdx, :)).^2, 2));
        distMatrix = min(distMatrix, distToNearestSelected);
        
        % Find the point with the maximum distance
        [~, maxIdx] = max(distMatrix);
        idx(i) = maxIdx;
        selectedIdx = maxIdx;
    end
end