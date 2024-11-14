% The MOMOS Model - Stripes
% Solving the PDE system with diffusion approximated by central finite differences
% and the chemotaxis term approximated by a second-order upwind scheme.
% Time-stepping is performed using the IMSP (IMplicit SymPlectic) first order method.

clear all
close all
clc

% PDE System
% u_t = Du * Lap(u) - div(beta * h(u) * grad(v)) + f(u,v)
% v_t = Dv * Lap(v) + g(u,v)

% Set parameters
Du = 1e-3;        % Diffusion coefficient for u
Dv = 1e-3;        % Diffusion coefficient for v
beta = 0.056;     % Chemotactic sensitivity
k1 = 0.4;         % Reaction rate constant
k2 = 0.6;         % Reaction rate constant
c = 1e-3;         % Constant in reaction term for v
q = 0.075;        % Reaction constant in f(u, v)

% Define reaction terms as functions
f = @(u, v) -k1*u-q*u.*abs(u)+k2*v;
g = @(u, v) k1*u-k2*v+c;

% Spatial domain
Lx = 25;          % Length of domain in x
Ly = 25;          % Length of domain in y
hx = 1;           % Grid spacing in x
hy = 1;           % Grid spacing in y
nx = round(Lx/hx) + 1;  % Number of grid points in x
ny = round(Ly/hy) + 1;  % Number of grid points in y

% Spatial grid
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[xx,yy] = meshgrid(x,y);

% Time domain
t0 = 0;
tf = 80000;
ht = 0.1;         % Time step size
t = t0:ht:tf;
nt = length(t);

% Spatial semi-discretization for diffusion with periodic boundary conditions
spatial_discretization_diffusion_periodic;

% Set initial conditions
ue = sqrt(c/q);
ve = (k1/k2)*sqrt(c/q) + c/k2;
U0 = ue * ones(nx, ny) + 1e-5 * rand(nx, ny);
V0 = ve * ones(nx, ny) + 1e-5 * rand(nx, ny);
% load('ic_momos_20.mat');  % Load custom initial conditions if available

% Initialize variables for solving the ODE system using IMSP scheme
u = U0(:);
v = V0(:);
u_sol(:,1) = u;
v_sol(:,1) = v;
mean_u(1) = mean(u);

% Precompute matrices for implicit scheme
In = eye(n);
Mu = inv(In-Du*ht*A);  % Matrix for implicit diffusion update in u
Mv = inv(In-Dv*ht*A);  % Matrix for implicit diffusion update in v

% Time-stepping loop
for k = 1 : nt-1
    % Update v with implicit step
    v = (v+ht*(k1*u+c))./(1+ht*k2);
    % Compute chemotaxis terms in x and y directions
    [Chemotaxis_x, Chemotaxis_y] = spatial_discretization_chemotaxis_upwind_second_order_periodic( ...
        reshape(u, ny, nx), reshape(v, ny, nx), hx, hy, nx, ny);
    % Update u with reaction and chemotaxis terms
    u = u + ht*f(u,v) - ht*beta*(Chemotaxis_x + Chemotaxis_y);
    u = Mu*u;  % Apply implicit diffusion update for u
    % Apply implicit diffusion update for v
    v = Mv*v;  
    % Flatten vectors for consistency
    u = u(:);
    v = v(:);
    % Store results
    u_sol(:,k+1) = u;
    v_sol(:,k+1) = v;
    mean_u(k+1) = mean(u);
end

% Plot final solution
figure
pcolor(xx,yy,reshape(u,ny,nx));
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
