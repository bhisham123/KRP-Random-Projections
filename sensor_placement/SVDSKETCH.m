function [U,inds,times] = SVDSKETCH(X,r,p,modes,method)
% Implelemtation of randomized HOSVD method based on SVDSKETCH suggested in the paper
% "Tensor-based flow reconstruction from optimally located sensor measurements"
% by M. Farazmand and A. K. Saibaba. Journal of Fluid Mechanics, 962:A27, 2023

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


% store dimensions and properties
inds = repmat({':'},length(modes),1);

t_rng = 0;
t_unf = 0;
t_mult = 0;
t_fact = 0;

% initialize factor matrices
U = cell(length(modes),1);

for n = modes
    tic;
    M = double(tenmat(X,n));
    % [Q,~,~] = svdsketch(M, 1.e-2); %'MaxSubspaceDimension'
    [Q,~,~] = svdsketch(M, 1.e-2,'MaxSubspaceDimension',r(n)+p);
    Q = Q(:,1:r(n)+p);
    inds{n} = subsetselection(Q,method);
    U{n} = Q/Q(inds{n},:);
    t_fact = t_fact + toc;
end
%timming summary
times = [t_mult, t_fact,t_rng,t_unf];
end