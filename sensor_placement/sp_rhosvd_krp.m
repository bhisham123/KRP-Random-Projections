function [U,inds,times] = sp_rhosvd_krp(X,r,p,modes,method)
% Structure-Preserving randomized HOSVD with KRP random matrix
%
% Inputs
%   X: original tensor (d modes)
%   r: target rank vector [r1,...,rd]
%   p: oversampling parameter
%   modes: Modes (dimensions) over which to compute factor matrices (subset of 1:d).
%   method: Column selection method used in 'passefficient' (e.g., 'pqr', 'rrqr' and 'deim')

% Outputs: 
%  U: Cell array of factor matrices for selected modes.
%  inds: Cell array of selected column indices for each mode (from passefficient).
%  times: 1×4 vector with timing information:
%                  [t_mttkrp, t_factor, t_random, t_unfold]


% Written by Bhisham Dev Verma, 2025

% store dimensions and properties
dims = size(X);
d = length(dims);

% Placeholder for selected indices (for column sampling)
inds = repmat({':'},length(modes),1);

% initialize factor matrices
U = cell(length(modes),1);

%initialize random matrices 
Phi = cell(d,1);   

% Timing variables
t_rng  = 0; % Time for random generation
t_unf  = 0; % Tensor unfolding time
t_fact = 0; % Time for computing factor matrices
t_mtt  = 0; % Time for MTTKRP computation


for n = modes
    %Generate random matrices
    tic;
    for j = [1:n-1,n+1:d]
        Phi{j} = randn(size(X,j),r(n)+p);
    end
    t_rng = t_rng + toc;
    
    % compute sketch
    tic;
    Y = mttkrp(X,Phi,n); t_mtt = t_mtt + toc;

    % compute factor matrix            
    tic;
    [U{n},inds{n}] = passefficient(Y,method); t_fact = t_fact + toc;
end
% timing summary
times = [t_mtt, t_fact,t_rng,t_unf];
end

