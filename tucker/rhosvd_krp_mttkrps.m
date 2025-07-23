function [T, times] = rhosvd_krp_mttkrps(X,r,p,modes)
% Randomized HOSVD with KRP random matrix based on memoization
%
% Inputs
%   X: original tensor (d modes)
%   r: target rank vector [r1,...,rd]
%   p: oversampling parameter
%   modes: order of processing the modes/dimensions (list)

% Outputs
%   T: Tucker tensor
%   time : 1×5 vector with timing information:
%                  [t_core, t_mtt, t_fact,t_rng, t_mat]


% Written by Bhisham Dev Verma, 2025

% store dimensions and properties
dims = size(X); %mode sizes
d = length(dims); %number of modes


t_rng  =  0;
t_mat = 0;
t_fact = 0;
t_mtt  =  0;
t_core = 0;


% Pre allocate memory
U = cell(length(modes),1);

%Generate random matrices
Phi = cell(1,d);
tic;
for n = 1:d
    Phi{n} = randn(size(X,n),r(n)+p);
end
t_rng = t_rng + toc;
    
% compute sketch
tic;
Y = mttkrps(X,Phi);
t_mtt = t_mtt + toc;
    
% compute factor matrix
tic;
for n = modes
    [U{n},~] = qr(Y{n},0);
end
t_fact = t_fact + toc;

%Compute core 
tic;
G = ttm(X,U,'t');
t_core = t_core + toc;

% return Tucker tensor
T = ttensor(G,U);
times = [t_core,t_mtt,t_fact,t_rng, t_mat];
end