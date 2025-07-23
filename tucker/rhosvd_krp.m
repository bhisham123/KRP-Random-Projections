function [T,times] = rhosvd_krp(X,r,p,modes)
% Randomized HOSVD with KRP random matrix
%
% Inputs
%   X: original tensor (d modes)
%   r: target rank vector [r1,...,rd]
%   p: oversampling parameter
%   modes: order of processing the modes (list)

% Outputs
%   T: Tucker tensor
%   time : 1×5 vector with timing information:
%                  [t_core, t_mtt, t_fact,t_rng, t_mat]


% Written by Bhisham Dev Verma, 2025


% mode size and number of modes
sz = size(X);
d = length(sz);
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

for n = modes
    tic;
    for j = [1:n-1,n+1:d]
        Phi{j} = randn(size(X,j),r(n)+p);
    end
    t_rng = t_rng + toc;
    
    tic;
    Y = mttkrp(X,Phi,n); t_mtt = t_mtt + toc;
                
    tic;
    [U{n},~,~] = qr(Y,0);  t_fact = t_fact + toc;
end
% Compute core tensor
tic;
G = ttm(X,U,'t');    t_core = t_core + toc;

% return Tucker tensor and time summary
T = ttensor(G,U);
times = [t_core, t_mtt, t_fact,t_rng, t_mat];
end

