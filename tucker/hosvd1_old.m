function [T,times] = hosvd1_old(X,r,modes)


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
    [Q,~,~] = svd(M,'econ');   % compute svd decomposition 
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
times = [t_core,t_mult,t_fact,0,t_mat];
end


