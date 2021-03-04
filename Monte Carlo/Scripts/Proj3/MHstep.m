function tnext = MHstep(t_prev, rho, fstat) %input conditional params
    
    %Draw X*
    [tstar, pdf] = rwp(t_prev, rho);
    prevPdf = prevpdf(t_prev, tstar, rho);
    
    %Trick to make sure that if the outer boundaries are passed the 
    %suggestion will be rejected
    if tstar(1) ~= t_prev(1) || tstar(end) ~= t_prev(end)
        prevPdf = 0;
    end
    
    %Calculation of alpha
    logExpression = log(fstat(tstar)*prevPdf)-log(fstat(t_prev))-log(pdf);
    alpha = min(1, exp(logExpression));
    U = rand;
    
    %Rejection step
    if U <= alpha
        tnext = tstar;
    else
        tnext = t_prev;
    end
    
end