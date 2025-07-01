function [interior_coords, K] = compute_gaussian_curvature(vertices, faces)
    % vertices - Nx3 matrix of vertex coordinates
    % faces - Mx3 matrix of triangle vertex indices
    
    num_vertices = size(vertices, 1);
    num_faces = size(faces, 1);
    
    % Initialize edge map
    edge_count = containers.Map('KeyType', 'char', 'ValueType', 'double');
    
    % Process each face to build edge count map
    for i = 1:num_faces
        % Get the vertices of the current face
        v1 = faces(i, 1);
        v2 = faces(i, 2);
        v3 = faces(i, 3);
        
        % Create edges (always sort vertices to ensure consistent key ordering)
        edges = [v1, v2; v2, v3; v3, v1];
        for j = 1:size(edges, 1)
            edge = sort(edges(j, :)); % Ensure consistent key ordering
            edge_key = sprintf('%d-%d', edge(1), edge(2));
            
            if isKey(edge_count, edge_key)
                edge_count(edge_key) = edge_count(edge_key) + 1;
            else
                edge_count(edge_key) = 1;
            end
        end
    end
    
    % Identify boundary edges (edges that appear only once)
    boundary_edges_keys = keys(edge_count);
    boundary_edges = [];
    for i = 1:length(boundary_edges_keys)
        if edge_count(boundary_edges_keys{i}) == 1
            edge = sscanf(boundary_edges_keys{i}, '%d-%d')';
            boundary_edges = [boundary_edges; edge];
        end
    end
    
    % Initialize boundary vertex list
    boundary_vertices = false(num_vertices, 1);
    
    % Mark boundary vertices based on boundary edges
    for i = 1:size(boundary_edges, 1)
        edge = boundary_edges(i, :);
        boundary_vertices(edge(1)) = true;
        boundary_vertices(edge(2)) = true;
    end

     % Find all triangles containing boundary vertices
    connected_faces = false(num_faces, 1);
    for i = 1:num_faces
        % Get the vertices of the current face
        v1 = faces(i, 1);
        v2 = faces(i, 2);
        v3 = faces(i, 3);
        
        % Check if any vertex of this face is a boundary vertex
        if boundary_vertices(v1) || boundary_vertices(v2) || boundary_vertices(v3)
            connected_faces(i) = true;
        end
    end
    
    % Find all vertices connected to boundary vertices through connected faces
    connected_vertices = false(num_vertices, 1);
    for i = 1:num_faces
        if connected_faces(i)
            % Get the vertices of this face
            v1 = faces(i, 1);
            v2 = faces(i, 2);
            v3 = faces(i, 3);
            
            % Mark vertices of this face as connected
            connected_vertices(v1) = true;
            connected_vertices(v2) = true;
            connected_vertices(v3) = true;
        end
    end
    
    % Update boundary vertices to include all connected vertices
    boundary_vertices = connected_vertices;
    
    % Identify interior vertices as those not marked as boundary vertices
    interior_vertices = ~boundary_vertices;
    
    % Initialize angle sums and counts
    angle_sums = zeros(num_vertices, 1);
    min_angles = inf(num_vertices, 1); % Track the minimum angle for each vertex
    max_angles = -inf(num_vertices, 1); % Track the maximum angle for each vertex
    
    % Compute angles for each triangle
    for i = 1:num_faces
        % Get the vertices of the current face
        v1 = faces(i, 1);
        v2 = faces(i, 2);
        v3 = faces(i, 3);
        
        % Get coordinates of vertices
        p1 = vertices(v1, :);
        p2 = vertices(v2, :);
        p3 = vertices(v3, :);
        
        % Compute vectors for each side of the triangle
        vec1 = p2 - p1;
        vec2 = p3 - p2;
        vec3 = p1 - p3;
        
        % Compute angles using dot product and vector magnitudes
        angle1 = acos(dot(vec1, -vec3) / (norm(vec1) * norm(vec3)));
        angle2 = acos(dot(vec2, -vec1) / (norm(vec2) * norm(vec1)));
        angle3 = acos(dot(vec3, -vec2) / (norm(vec3) * norm(vec2)));
        
        % Add angles to the sum for each vertex correctly
        angle_sums(v1) = angle_sums(v1) + angle1;
        angle_sums(v2) = angle_sums(v2) + angle2;
        angle_sums(v3) = angle_sums(v3) + angle3;
        
        % Track the minimum and maximum angle for each vertex
        angles = [angle1, angle2, angle3];
        min_angles(v1) = min(min_angles(v1), angles(1));
        min_angles(v2) = min(min_angles(v2), angles(2));
        min_angles(v3) = min(min_angles(v3), angles(3));
        max_angles(v1) = max(max_angles(v1), angles(1));
        max_angles(v2) = max(max_angles(v2), angles(2));
        max_angles(v3) = max(max_angles(v3), angles(3));
    end
    
    % Compute angle deficit for each vertex
    angle_deficit = 2 * pi - angle_sums;
    
    % Exclude vertices based on the following conditions:
    % 1. Angles smaller than 31° or larger than 89°
    % 2. Curvature greater than 2π/3
    exclude_vertices = (min_angles < deg2rad(33)) | (max_angles > deg2rad(87)) | (abs(angle_deficit)>0.01);
%     exclude_vertices_threshold = abs(angle_deficit) > (0.1);
%     exclude_vertices = exclude_vertices+exclude_vertices_threshold;
    interior_vertices = interior_vertices & ~exclude_vertices;
    
    % Assign curvature only to valid interior vertices
    K = zeros(num_vertices, 1);
    K(interior_vertices) = angle_deficit(interior_vertices);
    
    % Extract coordinates and curvatures for interior vertices
    interior_coords = vertices(interior_vertices, :);
    K = K(interior_vertices);
    
end
