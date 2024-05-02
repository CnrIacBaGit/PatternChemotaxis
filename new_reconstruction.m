function [error,sol_dmd] = new_reconstruction(dt,t0,T,W,L,X)
%% Ricostruzione Dinamica
 
 mu=diag(L);
 omega=log(mu)/dt;
 %y0=W\X(:,1);
 y0=pinv(W)*X(:,1);

t = t0:dt:T;

for j=1:length(t)
    sol_dmd(:,j)=y0.*exp(omega*t(j));
end

sol_dmd_full = real(W*sol_dmd);
%keyboard
%error = norm(sol_dmd_full-X,'fro')/norm(X,'fro');
error = max(max(abs(sol_dmd_full-X)))/max(max(abs(X)));

