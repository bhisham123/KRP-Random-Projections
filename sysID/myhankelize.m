function H = myhankelize(markov,p,q,no,ni)
    % Creates Hankel matrix from markov parameters
    % Inputs:
    %   markov: cell of markov parameters, each block no x ni,  
    %   p: number of block rows
    %   q: number of block columns

    H = zeros(p*no,q*ni);

    for i=1:p             
       B = cell2mat(markov(i:q+i-1));
       H((i-1)*no+1:i*no,:) = B;
    end
end

%H((i-1)*nr+1:i*nr,:) = markov(:,i*ns+1:(q+i)*ns);