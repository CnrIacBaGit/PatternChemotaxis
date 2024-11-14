% Spatial semi-discretization 
% Chemotacit term
% − div(u grad(v)) = −delta(u v) + div(v grad(u))

if  nx == ny
    n = nx^2;
    I = speye(nx);
    evec = ones(nx,1);
    % 1D Laplacian
    Lap = spdiags([evec,-2*evec,evec],[1,0,-1],nx,nx);
    % Homogeneous Neumann boundary conditions
    Lap(1,2) = 2;
    Lap(end,end-1) = 2;
    % 2D Laplacian in Kronecker form
    A = ((1/hx)^2)*kron(I,Lap) + ((1/hy)^2)*kron(Lap,I) ;
    % 1D nonlinear diffusion div(v grad(u))
    A1 = spdiags([evec,evec],[0,1],nx,nx);
    A2 = spdiags([evec,evec],[0,-1],nx,nx);
    B1 = spdiags([-evec,evec],[0,1],nx,nx);
    B2 = spdiags([evec,-evec],[0,-1],nx,nx);
    % Homogeneous Neumann boundary conditions
    A1(end,end-1) = 1;
    A2(1,2) = 1;
    B1(end,end-1) = 1;
    B2(1,2) = -1;
    A1 = 0.5*A1;
    A2 = 0.5*A2;
    % 2D nonlinear diffusion in Kronecker form
    A1x = (1/hx)*kron(I,A1);
    A1y = (1/hy)*kron(A1,I);
    A2x = (1/hx)*kron(I,A2);
    A2y = (1/hy)*kron(A2,I);
    B1x = (1/hx)*kron(I,B1);
    B1y = (1/hy)*kron(B1,I);
    B2x = (1/hx)*kron(I,B2);
    B2y = (1/hy)*kron(B2,I); 
else
    evecx = ones(nx,1);
    evecy = ones(ny,1);
    % 1D Laplacian along x and y
    Lapx = spdiags([evecx,-2*evecx,evecx],[1,0,-1],nx,nx);
    Lapy = spdiags([evecy,-2*evecy,evecy],[1,0,-1],ny,ny);
    % Homogeneous Neumann boundary conditions
    Lapx(1,2) = 2;
    Lapx(end,end-1) = 2;
    Lapy(1,2) = 2;
    Lapy(end,end-1) = 2;
    % 2D Laplacian in Kronecker form
    Inx = speye(nx);
    Iny = speye(ny);
    A = ((1/hx)^2)*kron(Iny, Lapx) + ((1/hy)^2)*kron(Lapy,Inx);
    % 1D nonlinear diffusion div(v grad(u)) along x and y
    A11 = spdiags([evecx,evecx],[0,1],nx,nx);
    A12 = spdiags([evecy,evecy],[0,1],ny,ny);
    A21 = spdiags([evecx,evecx],[0,-1],nx,nx);
    A22 = spdiags([evecy,evecy],[0,-1],ny,ny);
    B11 = spdiags([-evecx,evecx],[0,1],nx,nx);
    B12 = spdiags([-evecy,evecy],[0,1],ny,ny);
    B21 = spdiags([evecx,-evecx],[0,-1],nx,nx);
    B22 = spdiags([evecy,-evecy],[0,-1],ny,ny);
    % Homogeneous Neumann boundary conditions
    A11(end,end-1) = 1;
    A12(end,end-1) = 1;
    A21(1,2) = 1;
    A22(1,2) = 1;
    B11(end,end-1) = 1;
    B12(end,end-1) = 1;
    B21(1,2) = -1;
    B22(1,2) = -1;
    A11 = 0.5*A11;
    A12 = 0.5*A12;
    A21 = 0.5*A21;
    A22 = 0.5*A22;
    % 2D nonlinear diffusion in Kronecker form
    A1x = (1/hx)*kron(Iny,A11);
    A1y = (1/hy)*kron(A12,Inx);
    A2x = (1/hx)*kron(Iny,A21);
    A2y = (1/hy)*kron(A22,Inx);
    B1x = (1/hx)*kron(Iny,B11);
    B1y = (1/hy)*kron(B12,Inx);
    B2x = (1/hx)*kron(Iny,B21);
    B2y = (1/hy)*kron(B22,Inx); 
    n = nx*ny;
end

