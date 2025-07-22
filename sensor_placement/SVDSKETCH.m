function [U,inds,times] = SVDSKETCH(X,r,p,modes,method)

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
times = [t_mult, t_fact,t_rng,t_unf];
end