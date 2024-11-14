function [u_pdmd,v_pdmd,error,err_history,rango,store_rank] = pDMD(X,dt,int_1,t0,tol,tol2)
% pDMD function
dim = size(X,2);
n_int = size(X,2)/10; % max number of intervals (min 10 snapshots)
d_int = 1; % Increment for the intervals
success = [];
for j = int_1 : d_int : n_int
    j
    [j size(X,2)/j];
    dim_per = dim/j;
    sol_per = [];
    for i = 1 : j
        in = ceil((i-1)*dim_per+1);
        fine = ceil(i*dim_per);
        if i == j
            X_st = X(:,in:end);
        else
            X_st = X(:,in:fine);
        end
        target_rank(i) = rank(X_st);
        if target_rank(i)<=0
            target_rank(i) = 1;
        end
        err_qb_count = 1;
        while err_qb_count > tol
            tic;
            [W_qb, L_qb,A_qb,Q_qb] = rdmd_qb(X_st,target_rank(i));
            time_qb(i) =toc;
            [err_qb_count,sol_dmd_qb] = new_reconstruction(dt,t0,(size(X_st,2)-1)*dt,W_qb,L_qb,X_st);
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
            break;
        end
        err_qb(i) = err_qb_count;
        sol_per = [sol_per real(W_qb*sol_dmd_qb)];
        store_rank{j}(i) = target_rank(i);
        rango(j) = max(target_rank);
    end
    success = [success risultato];
    if success(end) >0
        error(j) = norm(sol_per-X,'fro')/norm(X,'fro');
        for k = 1:dim
            err_history{j}(k) = norm(X(:,k)-sol_per(:,k))/norm(X(:,k));
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
    clear time_qb
end
% pDMD solutions u and v
u_pdmd = sol_per(1:end/2,:);
v_pdmd = sol_per(end/2+1:end,:);
end

