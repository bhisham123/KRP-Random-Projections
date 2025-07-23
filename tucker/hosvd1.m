function [T,times] = hosvd1(X,r,modes)
% HOSVD method 

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
    tic;
     M = double(tenmat(X,n));
    t_mat = t_mat + toc;
    tic;
    [~,R] = qr(M','econ');  
    [Q,~,~] = svd(R',"econ");
    % Extract factor matrix by picking out leading eigenvectors 
    U{n} = Q(:,1:r(n));
    t_fact = t_fact+toc;
end
%Compute core 
tic;
 G = ttm(X,U,'t');
t_core = t_core + toc;

% return Tucker tensor
T = ttensor(G,U);
times = [t_core,t_mult,t_fact,0, t_mat];
end


