function [W L,A,Q] = rdmd_qb(X,perm)
% randomized DMD with QB decomposition
perm = perm;
[Q B] = rqb(X,perm);
B_in = B(:,1:end-1);
B_ou = B(:,2:end);
[U S V] = svds(B_in,perm);
A = U'*B_ou*V*pinv(S);
[W_tmp L] = eig(A);
W = Q*B_ou*V*pinv(S)*W_tmp;
A = B_ou*pinv(B_in);
end
