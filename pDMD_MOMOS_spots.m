% The application of pDMD for the MOMOS model
% spots

% Generate data
% We apply the pDMD on the time interval [t0,tbar] with tbar = 5000
MOMOS_spots_nld_pbcs_tf_5000

% Save the snapshots every 5 time steps
u_fom = u_sol(:,2:5:end);
v_fom = v_sol(:,2:5:end);
X = [u_fom; v_fom];
dt = 5*ht;
t0 = 0;
tf = tf - dt;

tol = 1e-1; % tol as in Algorithm 1 (step 8)
bar_tol = 1e-4; % tolerance Algorithm 1 (step 14)
N = 5; % Initial number of intervals
[u_pdmd,v_pdmd,error,err_history,rango,store_rank] = pDMD(X,dt,N,t0,tol,bar_tol);

% Plots pDMD
figure
pcolor(xx,yy,reshape(u_pdmd(:,end),ny,nx))
xlabel('x')
ylabel('y')
shading interp
colorbar
colormap('jet')
title('pDMD u at t = 5000')
set(gca,'FontSize',12,'FontWeight','B')

% Error
figure
semilogy(0:dt:(size(X,2)-1)*dt,err_history{end},'LineWidth',2)
xlabel('t')
title('Error over time')
set(gca,'FontSize',12,'FontWeight','B')
