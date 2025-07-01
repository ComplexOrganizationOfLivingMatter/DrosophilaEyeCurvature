function V_smoothed = LaplacianSmoothing(V, F, lambda, iterations)
    % V: Vertices (Nx3 matrix, N vertices in 3D)
    % F: Faces (Mx3 matrix, M triangular faces)
    % lambda: Smoothing factor (controls the step size)
    % iterations: Number of smoothing iterations
    
    N = size(V, 1);   % Number of vertices
    
    % Step 1: Build the adjacency matrix
    A = sparse(N, N); % Initialize sparse matrix
    for i = 1:size(F, 1)
        % For each triangle in F, connect the vertices
        A(F(i,1), F(i,2)) = 1;
        A(F(i,2), F(i,1)) = 1;
        A(F(i,2), F(i,3)) = 1;
        A(F(i,3), F(i,2)) = 1;
        A(F(i,3), F(i,1)) = 1;
        A(F(i,1), F(i,3)) = 1;
    end
    
    % Make A symmetric and binary (no repeated edges)
    A = (A + A') > 0;
    
    % Step 2: Compute the Laplacian matrix (L = D - A)
    D = diag(sum(A, 2)); % Degree matrix
    L = D - A;           % Laplacian matrix
    
    % Step 3: Iteratively smooth the mesh
    V_smoothed = V;  % Initialize the smoothed vertices
    for iter = 1:iterations
        % Perform one iteration of Laplacian smoothing
        V_smoothed = V_smoothed - lambda * (L * V_smoothed);
    end
end