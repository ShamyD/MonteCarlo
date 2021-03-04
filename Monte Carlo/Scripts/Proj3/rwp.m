function [candidate, pdf] = rwp(t_prev, rho)
%generates a random walk proposal step candidate for an MH-sampler 
%t_prev is a row vector - Not matrix
    len_t = length(t_prev);
    d1 = diag(-ones(len_t, 1));
    d2 = diag(ones(len_t, 1), -2);
    matr = d1(:,1:end-2) + d2(1:end-2,1:end-4);
    R = t_prev*matr*rho;
    candidate = t_prev(2:end-1) + (rand(1,len_t-2)-0.5)*2.*R;
    candidate = [t_prev(1),candidate, t_prev(end)];
    pdf = 1/prod(2*R);
end