function f = ftprod(t)
%returns the pdf value for a set of breakpoints
    d = length(t);
    matr = -diag(ones(1, d), 0)  + diag(ones(1,d-1), -1);
    matr = matr(:,1:end-1)
    f = prod(t*matr)
    
end

