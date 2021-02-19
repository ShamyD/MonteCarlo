function matr = multToMatr(mult)
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

