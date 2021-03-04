function pdf = prevpdf(t_prev, tstar, rho)
    len_t = length(tstar);
    d1 = diag(-ones(len_t, 1));
    d2 = diag(ones(len_t, 1), -2);
    matr = d1(:,1:end-2) + d2(1:end-2,1:end-4);
    R = tstar*matr*rho;
   
    
    interval = zeros(2, len_t-2);
    interval(1, :) = tstar(2:end-1) + R;
    interval(2, :) = tstar(2:end-1) - R;
%   
    logic1 = interval(2,:) <= t_prev(2:end-1);
    logic2 = interval(1,:) >= t_prev(2:end-1);
    logic = logic1.*logic2
    
    
    
    pdf = 1/prod(2*R)*prod(logic);
end

