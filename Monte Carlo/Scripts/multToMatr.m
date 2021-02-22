function matr = multToMatr(mult)
    % Generate matrix which after multiplication with walks matrix
    % resamples the walks matrix.
    % Takes in the multiplicity mult, which says how many should
    %survive.
    
    %If first walk should have 0 copies second row 2, and third row 1
    % copy matrix would start like this:
    
    %[ 0  1  0  0 ...
    %[ 0  1  0  0 ...
    %[ 0  0  1  0 ...
    %[ .
    %[ .
    
    N = length(mult);
    matr = zeros(N, N);
    ind = 1;
    for i = 1:N
        len = mult(i, 1);
%         if len ~= 0
        matr(ind:ind+len-1, i) = 1;
%         end
        
        ind = ind + len;
    end
    
    
end

