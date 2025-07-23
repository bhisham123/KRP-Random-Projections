function [T,times]  = sthosvd_old(X,r,modes)
% STHOSVD method based on svd

% Inputs
%   X: original tensor (d modes)
%   r: target rank vector [r1,...,rd]
%   modes: order of processing the modes/dimensions (list)

% Outputs
%   T: Tucker tensor
%   time : 1×5 vector with timing information:
%                  [t_core, t_mtt, t_fact,t_rng, t_mat]


% Written by Bhisham Dev Verma, 2025



sz = size(X); %mode sizes
d = length(sz); %number of modes

% Pre allocate memory for factor matrices
U = cell(length(modes),1);

t_mat = 0;
t_mult = 0;
t_fact = 0;
t_core = 0;

G = X;
for n = modes
    
    tic;
    M = double(tenmat(G,n));
    t_mat = t_mat + toc;

    tic;
    [Q,~,~] = svd(M,'econ');   % compute svd decomposition 
    % Extract factor matrix by picking out leading eigenvectors 
    U{n} = Q(:,1:r(n));
    t_fact = t_fact+toc;

    % shrink
    tic;
    G = ttm(G,U,n,'t');
    t_core = t_core + toc;
end
% return Tucker tensor and time summary
T = ttensor(G,U);
times = [t_core,t_mult,t_fact,0, t_mat];
end


