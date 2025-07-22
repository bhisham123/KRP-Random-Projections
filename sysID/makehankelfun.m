function HF = makehankelfun(H,p,q,no,ni,tsignal)
    N = p+q-1; 
    
    % tsignal transposes the matrix implicitly.
    if nargin < 6, tsignal = 'notransp'; end

    % Define FV and FVT up here
    
    FV = zeros(N,no*ni);
    FVT = zeros(N,no*ni); 
    for i = 1:N
        h = reshape(H{i},1,no*ni); 
        ht = reshape(H{i}',1,no*ni); 
    
        j = mod(p+i-1,N)+1; 
        jt = mod(q+i-1,N)+1; 
    
        FV(j,:) = h;
        FVT(jt,:) = ht; 
    end 
    
    FV = fft(FV); 
    FVT = fft(FVT);
    
    if strcmp(tsignal,'notransp')
        HF = @(X,tflag) doHmult(X,tflag,p,q,no,ni,FV,FVT);
    elseif strcmp(tsignal,'transp')
        HF = @(X,tflag) doHmult(X,tflag,q,p,ni,no,FVT,FV); 
    else
        error("Transposition not indicated properly!"); 
    end
      
            
end

function HX = doHmult(X,tflag,p,q,no,ni,FV,FVT)
        nk = size(X,2); 
        N  = p + q - 1; 
        
        if strcmp(tflag,'transp')
            temp = p;
            p = q;
            q = temp;

            temp2 = no;
            no = ni;
            ni = temp2;
        elseif not(strcmp(tflag,'notransp'))
            error("Must indicate whether matrix is transposed!"); 
        end
        
        FX = zeros(q, ni*nk);
        HX = zeros(p*no, nk);
        
            
        % FFT of the X data
        for idx = 0:(ni-1)
            FX(q:-1:1,(idx*nk+1):((idx+1)*nk)) = X((idx+1):ni:end,:);
        end
        FX = fft(FX,N); 

        % Compute the multiplication 
        for ixr = 0:(no-1)
            Z = zeros(N,nk); 
            for ixc = 0:(ni-1)
                
                % Check which FFT set to multiply with
                if strcmp(tflag,'transp')
                    fv = FVT( : , ixc*no + ixr + 1);
                else
                    fv = FV( : , ixc*no + ixr + 1);
                end
                
                %Increment the fft-domain matvec product
                Z  = Z + FX( : , (ixc*nk+1):(ixc*nk+nk)).*fv; 
            end
            
            % Compute ifft and assign to output
            Z = ifft(Z); 
            HX((ixr+1):no:end,:) = Z(1:p,:);   
        end
end