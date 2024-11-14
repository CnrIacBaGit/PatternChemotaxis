% Future prediction by pDMD

clear
clc
close all

prompt = "Choose model ";
% Use 'Mimura', 'MOMOS_spots' or 'MOMOS_stripes'
model = input(prompt)

switch model
    case 'Mimura'
        % Generate data in [t0,T] with T = 600
        Mimura_hexagons_nld_neumann;
        int_1 = 5; % Initial number of intervals for the pDMD
        d_int = 1; % Increment for the intervals
        % Save the snapshots every 5 time steps
        u_fom = u_sol(:,2:5:end);
        v_fom = v_sol(:,2:5:end);
        X = [u_fom; v_fom];
        dt = 5*ht;
        t0 = 0;
        tf = tf - dt;
        tol = 1e-1; % tol as in Algorithm 1 (step 8)
        bar_tol = 1e-5; % tol as in Algorithm 1 (step 14)
        tbar = 500;
        nt1 = length([t0:dt:tbar])-1; % Time steps for the reconstruction in [t0,tbar]
        % n_pred = 20000;
    case 'MOMOS_spots'
        % Generate data
        % Generate data in [t0,T] with T = 10000
        MOMOS_spots_nld_pbcs; % Periodic BCs
        int_1 = 20; 
        d_int = 1;
        % Save the snapshots every 5 time steps
        u_fom = u_sol(:,2:5:end);
        v_fom = v_sol(:,2:5:end);
        X = [u_fom; v_fom];
        dt = 5*ht;
        t0 = 0;
        tf = tf - dt;
        tol = 1e-1; % tol as in Algorithm 1 (step 8)
        bar_tol = 1e-4; % tolerance Algorithm 1 (step 14)
        tbar = 5000;
        nt1 = length([t0:dt:tbar])-1; % Time steps for the reconstruction in [t0,tbar]
        % n_pred = 100000;
    case 'MOMOS_stripes'
        % Generate data
        % Generate data in [t0,T] with T = 100000
        MOMOS_stripes_nld_pbcs % Periodic BCs
        int_1 = 5;
        d_int = 1;
        % Save the snapshots every 5 time steps
        u_fom = u_sol(:,2:5:end);
        v_fom = v_sol(:,2:5:end);
        X = [u_fom; v_fom];
        dt = 5*ht;
        t0 = 0;
        tf = tf - dt;
        tbar = 50000;
        tol = 1e-1; % tol as in Algorithm 1 (step 8)
        bar_tol = 1e-3; % tolerance Algorithm 1 (step 14)
        nt1 = length([t0:dt:tbar])-1; % Time steps for the reconstruction in [t0,tbar]
        % n_pred = 100000;
end

tt = [t0:dt:tf];
n_pred = length(tt)-nt1; % Time steps for the prediction 
% nt1 = length(tt)-n_pred;

% pDMD reconstruction of X in [t0,tbar] (time steps nt1)

[sol_per,N,W_qb,L_qb] = fun_pdmd_prediction(int_1,d_int,X,dt,t0,nt1,tol,bar_tol);

% Future prediction in (tbar,T]
tf = tt(end)-tt(nt1);
tf = tf-dt;
X_pred = X(:,nt1+1:end);
[err_qb_count,sol_dmd_qb] = new_reconstruction(dt,t0,tf,W_qb,L_qb,X_pred);
sol_dmd = [sol_per real(W_qb*sol_dmd_qb)];

errore_tf = norm(sol_dmd(:,end)-X(:,end),'fro')/norm(X(:,end),'fro');
errore = norm(sol_dmd-X,'fro')/norm(X,'fro');

nt_pred = size(X_pred,2);
sol_dmd_pred = W_qb*sol_dmd_qb;
for kk = 1 : nt_pred
    error_prediction(kk) = norm(sol_dmd_pred(:,kk)-X_pred(:,kk),'fro')/norm(X_pred(:,kk),'fro');
end