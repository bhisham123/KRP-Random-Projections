function [U,S,V] = blockstructuredsvd_krp(mat, r)
    %   BLOCKSTRUCTUREDSVD_KRP(mat, r) computes a rank-r 
    %   approximation of a structured matrix defined by block matrices mat.M.
    %   The method uses single view randomized SVD algorithm with KRP structure random matrices (Algorithm 4.1).
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
    
    
    % Range sketch
    ell =  2*r + 1;  
    % ell = r+20; 
    Vlst = randn(ell,size(mat.E{1},2)); 
    Wlst = randn(ell,size(mat.M{1},2));
    p = length(mat.E{1});
    Yt = zeros(ell,size(mat.E{1},1)*size(mat.M{1},1));
    
    m = size(mat.M{1},1);
    q = size(Vlst,2);
    for j = 1:p
        temp = Wlst*mat.M{j}';
        krp_temp = rowwise_kr(temp,Vlst(:,j:-1:1));
        Yt(:,1:j*m) =  Yt(:,1:j*m) + krp_temp;
    end

    for j = p+1:2*p-1
        temp = Wlst*mat.M{j}';
        krp_temp = rowwise_kr(temp,Vlst(:,q:-1:j-q+1));
        Yt(:,(j-p)*m+1:end) =  Yt(:,(j-p)*m+1:end) + krp_temp;
    end

    [Q,~] = qr(Yt', 0);
    
    
    % Corange sketch
    ellp = 2*ell + 1;  
    % ellp = ceil(1.5*ell);  
    Vlst = randn(ellp,size(mat.E{1},1)); 
    Wlst = randn(ellp,size(mat.M{1},1));
    Wt = zeros(ellp, size(mat.E{1},2)*size(mat.M{1},2));

    n = size(mat.M{1},2);
    q = size(Vlst,2);
    for j = 1:p
        temp1 = Wlst*mat.M{j};
        krp_temp1 = rowwise_kr(temp1,Vlst(:,j:-1:1));
        Wt(:,1:j*n) =  Wt(:,1:j*n) + krp_temp1;
    end
    for j = p+1:2*p-1
        temp1 = Wlst*mat.M{j};
        krp_temp1 = rowwise_kr(temp1,Vlst(:,q:-1:j-q+1));
        Wt(:,(j-p)*n+1:end) =  Wt(:,(j-p)*n+1:end) + krp_temp1;
    end
   
    W = Wt;
    
    % Compute low-rank approximation
    Psi = kr(Vlst', Wlst'); 
    X = (Psi'*Q)\W;
    [UX,S,V] = svd(X, 'econ');
    U = Q*UX;
    
    U = U(:,1:r);
    S = S(1:r, 1:r);
    V = V(:,1:r);
end