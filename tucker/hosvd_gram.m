function [T,times] = hosvd_gram(X,r,modes)
% HOSVD method based on gram approach

% Inputs
%   X: original tensor (d modes)
%   r: target rank vector [r1,...,rd]
%   modes: order of processing the modes/dimensions (list)

% Outputs
%   T: Tucker tensor
%   time : 1×5 vector with timing information:
%                  [t_core, t_mtt, t_fact,t_rng, t_mat]


% Written by Bhisham Dev Verma, 2025

% store dimensions and properties
sz = size(X);
d = length(sz);

% Pre allocate memory for factor matrices
U = cell(length(modes),1);

t_mat = 0;
t_mult = 0;
t_fact = 0;
t_core = 0;

for n = modes
    % Compute Gram matrix
    tic;
    Yk = double(tenmat(X,n));
    t_mat = t_mat + toc;
    
    tic;
    Z = Yk*Yk';
    t_mult = t_mult+toc;

    tic;
    % Compute eigenvalue decomposition
    [V,D] = eig(Z);
    [~,pi] = sort(diag(D),'descend');
    % Extract factor matrix by picking out leading eigenvectors of V
    U{n} = V(:,pi(1:r(n)));
    t_fact = t_fact+toc;
end
%Compute core 
tic;
G = ttm(X,U,'t');
t_core = t_core + toc;

% return Tucker tensor
T = ttensor(G,U);
times = [t_core,t_mult,t_fact,0,t_mat];
end


