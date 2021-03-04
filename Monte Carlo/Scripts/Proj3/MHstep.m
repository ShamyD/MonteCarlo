function [tnext, niNext,sum_rej] = MHstep(t_prev, rho, tau, lambdaPrev, ni, sum_rej, i) %input conditional params
    
    %Draw X*
    [tstar, pdf] = rwp(t_prev, rho);
    prevPdf = prevpdf(t_prev, tstar, rho);
    niStar = create_ni(tau,tstar);
    
    %Trick to make sure that if the outer boundaries are passed the 
    %suggestion will be rejected
    
    if tstar(1) ~= 1658 || tstar(end) ~= 1980
        prevPdf = 0;
    end
    
    %Calculation of alpha4
    logExpression = ftpost(tstar, lambdaPrev, niStar)+log(prevPdf)- ftpost(t_prev, lambdaPrev, ni)-log(pdf);
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