function [T,times] = rhosvd(X,r,p,modes)
% Randomized HOSVD 
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
dims = size(X);
d = length(dims);

% Pre allocate memory
U = cell(d,1);    

t_rng = 0;
t_mat = 0;
t_fact = 0;
t_core = 0;
t_mult = 0;



for n = modes
    % generate Gaussian matrix
    tic;
    Omega = randn(prod(dims(1:n-1))*prod(dims(n+1:end)), r(n)+p);
    t_rng = t_rng+toc;

    % compute sketch
    tic;
    M = double(tenmat(X,n));
    t_mat = t_mat + toc;

    tic;
    Y = M * Omega; t_mult = t_mult + toc;

    % compute factor matrix
    tic;
    [U{n},~] = qr(Y,0); t_fact = t_fact + toc;
end
%Compute core 
tic;
    G = ttm(X,U,'t');
t_core = t_core + toc;

% return Tucker tensor and time summary
T = ttensor(G,U);
times = [t_core,t_mult,t_fact,t_rng, t_mat];
end