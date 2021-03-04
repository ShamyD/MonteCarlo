function [tnext, niNext] = MHstep(t_prev, rho, tau, lambdaPrev, ni) %input conditional params
    
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
    
%     ftpost(tstar, lambdaPrev, niStar)
%     prevPdf
%     ftpost(t_prev, lambdaPrev, ni)
%     pdf
    '----------'
    logExpression = log(ftpost(tstar, lambdaPrev, niStar))+log(prevPdf)-log(ftpost(t_prev, lambdaPrev, ni))-log(pdf);
    e = exp(logExpression);
    alpha = min(1, e);
    U = rand;
    
    %Rejection step
    if U <= alpha
        tnext = tstar;
        niNext = niStar;
    else
        tnext = t_prev;
        niNext = ni;
    end
    
end