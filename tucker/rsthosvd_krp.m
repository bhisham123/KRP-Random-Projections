function [T,times] = rsthosvd_krp(X,r,p, modes)
% Randomized STHOSVD with KRP random matrix
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

sz = size(X); %mode sizes
d = length(sz); %number of modes
if length(r) == 1
    r = r*ones(d,1);
end

%Pre allocate memory
U = cell(1,d);     %Cell array for factor matirces
Phi = cell(1,d);   %Cell array for random matrices

t_rng = 0;
t_mat = 0;
t_mtt = 0;  
t_fact = 0;
t_core = 0;

G = X;
for n = modes
    tic;
    for j = [1:n-1,n+1:d]
        Phi{j} = randn(size(G,j),r(n)+p);
    end
    t_rng = t_rng + toc;
    
    tic;
    Y = mttkrp(G,Phi,n); t_mtt = t_mtt + toc;
                
    tic;
    [U{n},~] = qr(Y,0);  t_fact = t_fact + toc;
    
    % shrink
    tic;
    G = ttm(G,U,n,'t');    t_core = t_core + toc;
end
% return Tucker tensor and time summary
T = ttensor(G,U);
times = [t_core, t_mtt, t_fact,t_rng, t_mat];
end

