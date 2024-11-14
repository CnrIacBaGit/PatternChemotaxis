% The MOMOS model - spots
% Solve PDE system with diffusion and chemotaxis term
% periodic BCs
% The chemotaxis term is discretized as a nonlinear diffusion
% Looking for pattern formation by chemotaxis-driven instability
clear all
close all
clc

% PDE system
% u_t = Du*Lap(u) - div(beta h(u) grad(v))+ f(u,v)
% v_t = Dv*Lap(v) + g(u,v)

% Set parameters
Du = 0.6;
Dv = 0.6;
beta = 1.2;
k1 = 0.4;
k2 = 0.6;
c = 0.8;
q = 0.075;
% Reaction terms
p = 1;
f = @(u,v) -k1*u.^p-q*u.*u+k2*v;
g = @(u,v) k1*u-k2*v+c;
% Spatial domain
Lx = 25;
Ly = 25; 
hx = 0.5;
hy = 0.5;
nx = round(Lx/hx) + 1;
ny = round(Ly/hy) + 1;
% Spatial grid
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[xx,yy] = meshgrid(x,y);
% Time domain
t0 = 0;
tf = 10000;
ht = 0.01;
t = [t0:ht:tf];
nt = length(t);

% Spatial semi-discretization
spatial_discretization_periodic

% Set initial conditions
ue = sqrt(c/q);
ve = k1/k2*(c/q)^(p/2)+c/k2;
U0 = ue*ones(nx,ny) + 1e-5*rand(nx,ny);
V0 = ve*ones(nx,ny) + 1e-5*rand(nx,ny);

% Solve the resulting ODE system using the IMSP_IE scheme
u = U0(:);
v = V0(:);
u_sol(:,1) = u;
v_sol(:,1) = v;
mean_u(1) = mean(u);
In = eye(n);
Mu = inv(In-Du*ht*A);
Mv = inv((1+k2*ht)*In-Dv*ht*A);
for k = 1 : nt-1
    k*ht
    v = Mv*(v+ht*(k1*u+c));
    u = u + ht*f(u,v) - ht*beta*((A1x*u).*(B1x*v)-(A2x*u).*(B2x*v)+(A1y*u).*(B1y*v)-(A2y*u).*(B2y*v));
    u = Mu*u;
    u_sol(:,k+1) = u;
    v_sol(:,k+1) = v;
    mean_u(k+1) = mean(u);
end

% Plots
figure
pcolor(xx,yy,reshape(u,ny,nx))
xlabel('x')
ylabel('y')
shading interp
colorbar
colormap('jet')
title('u at T = 10000')
set(gca,'FontSize',12,'FontWeight','B')

figure
plot(t,mean_u)
title('<u(t)>')
xlabel('t')
set(gca,'FontSize',12,'FontWeight','B')

for k = 1 : nt-1
    error(k) = norm(u_sol(:,k+1)-u_sol(:,k));
end

figure
semilogy(t(2:end),error)
xlabel('t')
set(gca,'FontSize',12,'FontWeight','B')
title('$\|u_{n+1}-u_n\|_2$','Interpreter','latex','FontSize',16)
