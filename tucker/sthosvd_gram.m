function [T,times]  = sthosvd_gram(X,r,modes)

sz = size(X);
d = length(sz);

% Pre allocate memory for factor matrices
U = cell(length(modes),1);

t_mat = 0;
t_mult = 0;
t_fact = 0;
t_core = 0;

G = X;
for n = modes
    tic;
    % Compute Gram matrix
    Yk = double(tenmat(G,n));
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

    % shrink
    tic;
    G = ttm(G,U,n,'t');
    t_core = t_core + toc;
end
T = ttensor(G,U);
times = [t_core,t_mult,t_fact,0, t_mat];
end


