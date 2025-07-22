function [A,B,C,t] = impulse_era(markov,p,q,no,ni,r,method)

    %rsvd parameters
    rho = 20;
    maxiter = 0;


    switch method
        case 'full'
            tic
            H = hankelize(markov,p,q,no,ni);
           [U,S,V] = svd(H,0);
            t = toc;

        case 'svds'
            tic
            H = hankelize(markov,p,q,no,ni);
            [U,S,V] = svds(H,r);
            t = toc;

        case 'svdshankel'
            tic
            H = cell2mat(markov);
            T = zeros(no*ni,p+q-1);   
            for i = 1:ni
                T(i:ni:end,:) = H(1:no,i:ni:end);
            end
            [U,S,V] = svds(@(x,tflag) svds_blockhank(x,tflag,T,no,ni,p,q),[p*no,q*ni],r);
            t = toc;
        case 'randsvd'
            t5 = tic;
            H = myhankelize(markov,p,q,no,ni);
            % size(H)
            [U,S,V] = rsvd(H,r,rho,maxiter);
            t = toc(t5);

        case 'singleview'
            Hf = makehankelfun(markov,p,q,no,ni,'notransp');
            Omega = randn(p*no,2*r+1);
            Psi = randn( q*ni, 4*r+2);
            Y = Hf(Omega, 'transp'); [Q,~] = qr(Y,0);
            W = Hf(Psi, 'notransp');
            X = (Psi'*Q)\W';
            [UX,S,V] = svd(X, 0);
            U = Q*UX;
    
            U = V(:,1:r);
            S = S(1:r, 1:r);
            V = U(:,1:r);
        
        case 'randkrp'
            Hf = makehankelfun(markov,p,q,no,ni,'notransp');
           
            Omega = kr(randn(q,r+p), randn(ni,r+p));
            Y = Hf(Omega, 'notransp'); [Q,~] = qr(Y,0);
            W = Hf(Q, 'transp');
            [UW,S,V] = svd(W','econ');
            U = Q*UW;
            
            
        case 'randsvdhankel'
            tic
            %Hf = blockhankel(cell2mat(markov),p,q,no,ni,0);
            Hf = makehankelfun(markov,p,q,no,ni,'notransp');
            [U,S,V] = rsvdFun(Hf,q*ni,r,rho,maxiter);
            t = toc;

        case 'tera'

            %tera parameters
            epsilon = 1e-2;

            N = size(markov,2);
            [markov_proj, W1, W2,~] = tera(markov,epsilon);  %W1 is no x r1, W2 is ni x r2       
            H = cell2mat(markov_proj(1:N));
            [U,S,V] = svd(H,'econ');

            no = size(W1,2); %r1
            ni = size(W2,2); %r2
            t=0;

        case 'randtera'
            N = size(markov,2);
            epsilon = 1e-2;
            [markov_proj, W1, W2,t_newm,s1,s2] = tera(markov,epsilon);  %W1 is no x r1, W2 is ni x r2
            no = size(W1,2); %r1
            ni = size(W2,2); %r2
            
            [no*p,ni*q]

            tic; H = cell2mat(markov_proj(1:N));
            Hf = blockhankel(H,p,q,no,ni);
            [U,S,V] = rsvd(Hf,r,rho,maxiter);
            t_r = toc;
            t = t_newm+t_r;
    end

    [A,B,C] = era(U,S,V,r,ni,no);

    if strcmp(method,'randtera')|| strcmp(method,'tera')
        B = B*W2';
        C = W1*C;
    end

end
