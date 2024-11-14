function [sol_per,N,W_qb,L_qb] = fun_pdmd_prediction(int_1,d_int,X,dt,t0,nt1,tol,tol2)
% Snapshot matrix X corresponding to the overall time interval t = [t0:dt:tf]
% This function computes the pDMD in the time interval [t0,tbar] 
% that corresponds to t(1:nt1)
X_train = X(:,1:nt1);
n_int = size(X_train,2)/10; % maximum number of intervals (min 10 snapshots)
dim = size(X_train,2);
success =[];
sol_per = [];
error = 0;
for j = int_1 : d_int : n_int
    [j size(X_train,2)/j]
    dim_per = dim/j;
    sol_per = [];
    for i = 1:j
        in = ceil((i-1)*dim_per+1);
        fine = ceil(i*dim_per);
        if i == j % The last time subinterval used for the prediction
            X_st = X_train(:,in:end);
            tf = (size(X_st,2)-1)*dt;
        else
            X_st = X_train(:,in:fine);
            tf = (size(X_st,2)-1)*dt;
        end
        target_rank(i) = rank(X_st);
        if target_rank(i)<=0
            target_rank(i) = 1;
        end
        err_qb_count = 1;        
        while err_qb_count > tol
            tic;
            [W_qb,L_qb,A_qb,Q_qb] = rdmd_qb(X_st,target_rank(i));
            time_qb(i) =toc;
            [err_qb_count,sol_dmd_qb] = new_reconstruction(dt,t0,tf,W_qb,L_qb,X_st);
            if isnan(err_qb_count) == 1
                err_qb_count = 1;
            end
            target_rank(i) = target_rank(i) - 1;
            if target_rank(i)<=0
                risultato = -1;
                break;
            else
                risultato = j;
            end
        end
        if risultato <0
            disp('pDMD does not converge')
            break;
        end
        err_qb(i) = err_qb_count;
        sol_per = [sol_per real(W_qb*sol_dmd_qb)];
        store_rank{j}(i) = target_rank(i); 
        rango(j) = max(target_rank);
    end
    success = [success risultato];
    if success(end) >0
        error(j) = norm(sol_per-X_train,'fro')/norm(X_train,'fro');
        for k = 1:dim
            err_history{j}(k) = norm(X_train(:,k)-sol_per(:,k))/norm(X_train(:,k));
        end
        tempo = sum(time_qb);
        time(j) = tempo;
        if error(end)< tol2
            break;
        end
    else
        time(j) = 0;
        error(j) = -100;
    end
    error(end);
    clear time_qb
end

N = length(error);