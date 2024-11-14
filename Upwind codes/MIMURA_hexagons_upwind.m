% The Mimura-Tsujikawa Model - Hexagons
% Solving the PDE system with diffusion approximated by central finite differences
% and the chemotaxis term approximated by a second-order upwind scheme.
% Time-stepping is performed using the Symplectic Euler method.

clear all
close all
clc

% PDE System
% u_t = Du * Lap(u) - div(beta * h(u) * grad(v)) + f(u, v)
% v_t = Dv * Lap(v) + g(u, v)

% Set parameters
Du = 0.0625;     % Diffusion coefficient for u
Dv = 1;          % Diffusion coefficient for v
beta = 17;       % Chemotactic sensitivity parameter
q = 7;           % Reaction constant in f(u, v)
k1 = 1;          % Reaction rate constant
k2 = 32;         % Reaction rate constant

% Define reaction terms as functions
f = @(u,v) q.*u.*(1-u);
g = @(u,v) k1.*u - k2.*v;

% Spatial domain
Lx = 3;          % Length of domain in x
Ly = 3;          % Length of domain in y
nx = 50;         % Number of grid points in x
ny = nx;         % Number of grid points in y (square grid)
hx = Lx/(nx-1);  % Grid spacing in x
hy = Ly/(ny-1);  % Grid spacing in y

% Spatial grid
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[xx, yy] = meshgrid(x,y);

% Time domain
t0 = 0;
tf = 600;        % Final time
ht = 1e-3;       % Time step size
t = t0:ht:tf;
nt = length(t);

% Spatial semi-discretization for diffusion with zero Neumann boundary conditions
spatial_discretization_diffusion_neumann;

% Set initial conditions
ue = 1;
ve = k1/k2;
U0 = ue*ones(nx,ny) + 1e-5*rand(nx,ny);
V0 = ve*ones(nx,ny) + 1e-5*rand(nx,ny);
% load("ic_kuto.mat");  % Load predefined initial conditions if available

% Initialize variables for solving the ODE system using the Symplectic Euler scheme
u = U0(:);
v = V0(:);
u_sol(:,1) = u;
v_sol(:,1) = v;
mean_u(1) = mean(u);

% Precompute matrices for implicit scheme
In = eye(n);
Mv = inv((1+ht*k2)*In-Dv*ht*A);  % Matrix for implicit diffusion update in v

% Time-stepping loop
for k = 1 : nt-1
    % Update v with implicit step
    v = Mv*(v+ht*k1.*u);
    % Compute chemotaxis terms in x and y directions
    [Chemotaxis_x,Chemotaxis_y] = spatial_discretization_chemotaxis_upwind_second_order_neumann( ...
        reshape(u,ny,nx),reshape(v,ny,nx),hx,hy,nx,ny);
    % Update u with reaction, diffusion, and chemotaxis terms
    u = u + ht*f(u,v) + ht*Du*A*u - ht*beta*(Chemotaxis_x + Chemotaxis_y);
    % Store results
    u_sol(:,k+1) = u;
    v_sol(:,k+1) = v;
    mean_u(k+1) = mean(u);
end

% Plot final solution for u
figure
pcolor(xx,yy,reshape(u,ny,nx));
xlabel('x');
ylabel('y');
shading interp
colorbar
colormap('jet')


% Plot mean u over time
figure
plot(t,mean_u)
title('<u(t)>')
xlabel('t')

% Optional: Generate animation
% h = figure;
% h.Visible = 'off';
% count = 0;
% for k = 1:100:size(u_sol, 2)
%     count = count + 1;
%     pcolor(xx, yy, reshape(u_sol(:, k), nx, ny));
%     shading interp;
%     colormap('jet');
%     drawnow;
%     M(count) = getframe;
% end
% h.Visible = 'on';
% movie(M);
