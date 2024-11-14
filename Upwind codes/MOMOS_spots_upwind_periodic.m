% The MOMOS Model - Spots
% Solving the PDE system with diffusion approximated by central finite differences
% and the chemotaxis term approximated by a second-order upwind scheme.
% Time-stepping is performed using the IMSP_IE method.

clear all
close all
clc

% PDE System
% u_t = Du * Lap(u) - div(beta * h(u) * grad(v)) + f(u, v)
% v_t = Dv * Lap(v) + g(u, v)

% Set parameters
Du = 0.6;        % Diffusion coefficient for u
Dv = 0.6;        % Diffusion coefficient for v
beta = 1.25;     % Chemotactic sensitivity parameter
k1 = 0.4;        % Reaction rate constant
k2 = 0.6;        % Reaction rate constant
c = 0.8;         % Constant in reaction term for v
q = 0.075;       % Reaction constant in f(u, v)

% Define reaction terms as functions
f = @(u, v) -k1*u - q*u.*abs(u) + k2*v;
g = @(u, v) k1*u - k2*v + c;

% Spatial domain
Lx = 25;         % Length of domain in x
Ly = 25;         % Length of domain in y
hx = 0.5;        % Grid spacing in x
hy = 0.5;        % Grid spacing in y
nx = round(Lx/hx) + 1;  % Number of grid points in x
ny = round(Ly/hy) + 1;  % Number of grid points in y

% Spatial grid
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[xx,yy] = meshgrid(x,y);

% Time domain
t0 = 0;
tf = 10000;     % Final time
ht = 0.01;      % Time step size
t = t0:ht:tf;
nt = length(t);

% Spatial semi-discretization for diffusion with zero Neumann boundary conditions
spatial_discretization_diffusion_periodic;

% Set initial conditions
ue = sqrt(c/q);
ve = (k1/k2)*sqrt(c/q) + c/k2;
U0 = ue*ones(nx,ny) + 1e-5*rand(nx,ny);
V0 = ve*ones(nx,ny) + 1e-5*rand(nx,ny);

% Solve the resulting ODE system using the IMSP_IE scheme
u = U0(:);
v = V0(:);
u_sol(:,1) = u;
mean_u(1) = mean(u);

% Precompute matrices for implicit scheme
In = eye(n);
Mu = inv(In-Du*ht*A);  % Matrix for implicit diffusion update in u
Mv = inv((1+k2*ht)*In-Dv*ht*A);  % Matrix for implicit diffusion update in v

% Time-stepping loop
for k = 1 : nt-1
    % Update v with implicit step
    v = Mv*(v+ht*(k1*u+c));
    % Compute chemotaxis terms in x and y directions
    [Chemotaxis_x, Chemotaxis_y] =  spatial_discretization_chemotaxis_upwind_second_order_periodic( ...
        reshape(u, ny, nx), reshape(v, ny, nx), hx, hy, nx, ny);
    % Update u with reaction and chemotaxis terms
    u = u + ht*f(u,v) - ht*beta*(Chemotaxis_x + Chemotaxis_y);
    u = Mu*u;  % Apply implicit diffusion update for u
    % Flatten vectors for consistency
    u = u(:);
    u_sol(:, k+1) = u;
    mean_u(k+1) = mean(u);
end

% Plot final solution for u
figure
pcolor(xx,yy,reshape(u,ny,nx))
xlabel('x')
ylabel('y')
shading interp
colorbar
colormap('jet')

% Plot mean u over time
figure
plot(t,mean_u)
title('<u(t)>')
xlabel('t')

% Optional: Generate animation (uncomment to visualize)
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
