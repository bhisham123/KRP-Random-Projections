function X = rowwise_kr(U, varargin)
%ROWWISE_KR Row-wise Kronecker product using bsxfun
%   X = rowwise_kr(A, B, ...) or rowwise_kr({A, B, ...})
%   Computes row-wise Kronecker product: kron(A(i,:), B(i,:), ...)
%   for each row i.

if ~iscell(U), U = [{U} varargin];
else U = [U varargin];
end

% Number of rows
R = size(U{1}, 1);

% Check all matrices have same number of rows
if any(cellfun(@(x) size(x,1), U) ~= R)
    error('All input matrices must have the same number of rows.');
end

X = U{1};  % Start with the first matrix

for n = 2:length(U)
    A = X;         % [R × I]
    B = U{n};      % [R × J]
    I = size(A,2);
    J = size(B,2);

    % Reshape A to [R, I, 1], B to [R, 1, J]
    A_exp = reshape(A, [R, I, 1]);
    B_exp = reshape(B, [R, 1, J]);

    % bsxfun applies row-wise outer product
    C = bsxfun(@times, A_exp, B_exp);   % Result: [R, I, J]

    % Reshape to [R, I*J]
    X = reshape(C, [R, I*J]);
end
