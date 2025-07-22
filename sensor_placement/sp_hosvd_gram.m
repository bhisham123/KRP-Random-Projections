function [U,inds,times] = sp_hosvd_gram(X,r,modes,method)
% Structure-Preserving HOSVD using Gram approach
%
% Inputs
%   X: original tensor (d modes)
%   r: target rank vector [r1,...,rd]
%   modes: Modes (dimensions) over which to compute factor matrices (subset of 1:d).
%   method: Column selection method used in 'subsetselection' (e.g., 'pqr', 'rrqr' and 'deim')

% Outputs: 
%  U: Cell array of factor matrices for selected modes.
%  inds: Cell array of selected column indices for each mode (from passefficient).
%  times: 1×4 vector with timing information:
%                  [t_mttkrp, t_factor, t_random, t_unfold]


% Written by Bhisham Dev Verma, 2025


% Placeholder for selected indices (for column sampling)
inds = repmat({':'},length(modes),1);

% Timing variables
t_rng  = 0;  % Time for random generation
t_unf  = 0;  % Tensor unfolding time
t_mult = 0;  % Time for Gram computation
t_fact = 0;  % Time for computing factor matrices

% initialize factor matrices
U = cell(length(modes),1);

for n = modes
    %Tensor Unfolding
    tic;
    Yk = double(tenmat(X,n));
    t_unf = t_unf + toc;

    % Compute Gram matrix
    tic;
    Z = Yk*Yk';
    t_mult = t_mult+toc;

    tic;
    % Compute eigenvalue decomposition
    [V,D] = eig(Z);
    [eigvec,pi] = sort(diag(D),'descend');
    % Extract factor matrix by picking out leading eigenvectors of V
    Q = V(:,pi(1:r(n)));
    inds{n} = subsetselection(Q,method);
    U{n} = Q/Q(inds{n},:);
    t_fact = t_fact + toc;
end
%timming sumary
times = [t_mult, t_fact,t_rng,t_unf];
end
