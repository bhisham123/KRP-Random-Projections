function [U,S,V] = rsvd(A,k,p,maxiter)
% A: funmat representing multiplication of any matrix and transpose
% k: target rank
% p: oversampling parameter

n = size(A,2);
Omega = randn(n,k+p);

%% Initial iteration
Y = A*Omega;
[Q,~] = qr(Y,0);

%% Subsequent iterations
    for j = 1:maxiter
        Y = A'*Q;
        [Q,~]=qr(Y,0);
        Y = A*Q;
        [Q,~]=qr(Y,0);
    end
    
%% Low rank approximation
    B = A'*Q;
    B = B';
    
    [U,S,V] = svd(B,'econ');
    U = Q*U;
    
 %% Compress
    U = U(:,1:k);
    V = V(:,1:k);
    S = S(1:k,1:k);
end