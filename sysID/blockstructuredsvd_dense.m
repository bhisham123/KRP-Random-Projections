function [U,S,V] = blockstructuredsvd_dense(mat, r)
    %   BLOCKSTRUCTUREDSVD_DENSE(mat, r) computes a rank-r 
    %   approximation of a structured matrix defined by block matrices mat.M.
    %   The method uses single view randomized SVD algorithm with dense Gaussian matrices.
    %
    %   Inputs:
    %     mat - Struct with fields:
    %           mat.M: cell array of 2p-1 dense blocks (m x n each)
    %           mat.E: used for inferring block sizes
    %     r   - Target rank for the approximation
    %
    %   Outputs:
    %     U - Approximate left singular vectors (size m*p x r)
    %     S - Diagonal matrix of top-r singular values (r x r)
    %     V - Approximate right singular vectors (size n*q x r)

    % Written by Bhisham Dev Verma, 2025

    
    % Range sketch
    % ell  =  2*r + 1; 
    ell = r+20;
    [p,q] = size(mat.E{1});
    [m,n] = size(mat.M{1});
    Omega = randn(n*q,ell); 
    
    Y = zeros(size(mat.E{1},1)*size(mat.M{1},1),ell);
    for k = 1:ell
        Omega_reshaped = reshape(Omega(:,k),n,q);
        for j = 1:p
            Omega_reshaped_trunct = Omega_reshaped(:,j:-1:1);
            Z = mat.M{j} * Omega_reshaped_trunct;
            Y(1:j*m,k) = Y(1:j*m,k)+reshape(Z,[],1);
        end

        for j = p+1:2*p-1
            Omega_reshaped_trunct = Omega_reshaped(:,q:-1:j-q+1);
            Z = mat.M{j} * Omega_reshaped_trunct;
            Y((j-p)*m+1:end,k) = Y((j-p)*m+1:end,k)+reshape(Z,[],1);
        end
    end


    [Q,~] = qr(Y, 0);
    
    
    % Corange sketch
    % ellp = 2*ell + 1; 
    ellp = ceil(1.5*ell); 
    Psi = randn(m*p,ellp); 
    W = zeros(n*q,ellp);
    for k = 1:ellp
        Psi_reshaped = reshape(Psi(:,k),m,p);
        for j = 1:p
            Psi_reshaped_trunct = Psi_reshaped(:,j:-1:1);
            Z = mat.M{j}' * Psi_reshaped_trunct;
            W(1:j*n,k) = W(1:j*n,k)+reshape(Z,[],1);
        end

        for j = p+1:2*p-1
            Psi_reshaped_trunct = Psi_reshaped(:,q:-1:j-q+1);
            Z = mat.M{j}' * Psi_reshaped_trunct;
            W((j-p)*n+1:end,k) = W((j-p)*n+1:end,k)+reshape(Z,[],1);
        end
    end
    W = W';
    
    % Compute low-rank approximation
    X = (Psi'*Q)\W;
    [UX,S,V] = svd(X, 'econ');
    U = Q*UX;
    
    U = U(:,1:r);
    S = S(1:r, 1:r);
    V = V(:,1:r);
end