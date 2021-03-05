function [tnext, niNext,sum_rej] = MHstep(t_prev, rho, tau, lambdaPrev, ni, sum_rej, i) %input conditional params
    
    %Draw X*
    [tstar, pdf, rev_pdf] = rwp(t_prev, rho);
    niStar = create_ni(tau,tstar);
    
    %Calculation of alpha
    logExpression = ftpost(tstar, lambdaPrev, niStar)+log(rev_pdf)- ftpost(t_prev, lambdaPrev, ni)-log(pdf);
    e = exp(logExpression);
    alpha = min(1, e);
    U = rand;
    
    %Rejection step
    if U <= alpha
        tnext = tstar;
        niNext = niStar;
        i;
    else
        tnext = t_prev;
        niNext = ni;
        sum_rej = sum_rej + 1;
    end
    
end