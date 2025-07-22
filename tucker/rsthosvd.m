function [T,times] = rsthosvd(X,r,p,modes)


% store dimensions and properties
dims = size(X);
d = length(dims);



%Pre allocate memory
U = cell(d,1);    

t_rng = 0;
t_mat = 0;
t_fact = 0;
t_core = 0;
t_mult = 0;


G = X;
for n = modes
    % generate Gaussian matrix
    tic;
    Omega = randn(prod(r(1:n-1))*prod(dims(n+1:end)), r(n)+p);
    t_rng = t_rng+toc;

    % compute sketch
    tic;
    M = double(tenmat(G,n));
    t_mat = t_mat + toc;

    tic;
    Y = M * Omega; t_mult = t_mult + toc;

    % compute factor matrix
    tic;
    [U{n},~] = qr(Y,0); t_fact = t_fact + toc;
    
    % shrink
    tic;
    G = ttm(G,U,n,'t');  t_core = t_core + toc;
end
% return Tucker tensor
T = ttensor(G,U);
times = [t_core,t_mult,t_fact,t_rng, t_mat];
end