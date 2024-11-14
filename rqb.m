function [Q B] = rqb(X,perm)
% QB decomposition
[n,m] = size(X);
Omega = randn(m,perm);
Y = X*Omega;
for i =1:2
    [Q, ~] = qr(Y,0);
    [Z, ~] = qr(X'*Q,0);
    Y = X*Z;
end
[Q, R] = qr(Y,0);
B = Q'*X;

