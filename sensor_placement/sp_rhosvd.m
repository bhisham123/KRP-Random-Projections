function [U,inds,times] = sp_rhosvd(X,r,p,modes,method)
% Structure-Preserving randomized HOSVD with dense Gaussian random matrix
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


% store dimensions 
dims = size(X);

% Placeholder for selected indices (for column sampling)
inds = repmat({':'},length(modes),1);

% Timing variables
t_rng  = 0;  % Time for random generation
t_unf  = 0;  % Tensor unfolding time
t_mult = 0;  % Time for unfolded tensor times matrix multiplication
t_fact = 0;  % Time for computing factor matrices

% initialize factor matrices
U = cell(length(modes),1);

for n = modes   
    % Generate random matrices
    tic;
    Phi = randn(prod(dims(1:n-1))*prod(dims(n+1:end)),r(n)+p);
    t_rng = t_rng + toc;
   
    % Tensor unfolding
    tic;
    M = double(tenmat(X,n));
    t_unf = t_unf + toc;
    
    % compute sketch
    tic;
    Y = M * Phi; t_mult = t_mult + toc;
    
    % compute factor matrix            
    tic;
    [U{n},inds{n}] = passefficient(Y,method); t_fact = t_fact+ toc;    
end
%timming summary
times = [t_mult,t_fact,t_rng,t_unf];
end