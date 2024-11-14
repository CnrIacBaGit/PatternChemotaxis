% Spatial semi-discretization for diffusion with zero Neumann boundary conditions

if nx == ny
    % Square grid case (nx = ny)
    n = nx^2;  % Total number of grid points
    I = speye(nx);  % Identity matrix of size nx
    evec = ones(nx, 1);  % Vector of ones for constructing Laplacian
    % 1D Laplacian matrix
    Lap = spdiags([evec, -2*evec, evec], [1, 0, -1], nx, nx);
    % Apply homogeneous Neumann boundary conditions
    Lap(1, 2) = 2;         % Adjust boundary point at the beginning
    Lap(end, end-1) = 2;   % Adjust boundary point at the end
    % 2D Laplacian in Kronecker form
    A = ((1/hx)^2) * kron(I, Lap) + ((1/hy)^2) * kron(Lap, I);
else
    % Non-square grid case (nx ≠ ny)
    evecx = ones(nx, 1);  % Vector of ones for x-direction Laplacian
    evecy = ones(ny, 1);  % Vector of ones for y-direction Laplacian
    % 1D Laplacian matrices along x and y directions
    Lapx = spdiags([evecx, -2*evecx, evecx], [1, 0, -1], nx, nx);
    Lapy = spdiags([evecy, -2*evecy, evecy], [1, 0, -1], ny, ny);
    % Apply homogeneous Neumann boundary conditions
    Lapx(1,2) = 2;            % Adjust boundary point at the beginning for Lapx
    Lapx(end,end-1) = 2;      % Adjust boundary point at the end for Lapx
    Lapy(1,2) = 2;            % Adjust boundary point at the beginning for Lapy
    Lapy(end,end-1) = 2;      % Adjust boundary point at the end for Lapy
    % 2D Laplacian in Kronecker form
    Inx = speye(nx);  % Identity matrix of size nx for x-direction
    Iny = speye(ny);  % Identity matrix of size ny for y-direction
    A = ((1/hx)^2)*kron(Iny,Lapx) + ((1/hy)^2)*kron(Lapy, Inx); 
    n = nx*ny;  % Total number of grid points for non-square grid
end