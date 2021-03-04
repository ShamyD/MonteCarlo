function tnext = MHstep(t_matr, rho, fstat, prevPdf)
    
    %Draw X*
    t_prev = t_matr(end, :)
    [tstar, pdf] = rwp(t_prev, rho)
    
    %Does this step have to be logged?
    alpha = min(1, fstat(tstar)*prevPdf/(fstat(t_prev)*pdf));
    U = rand;
    
    if U < alpha
        tnext = tstar;
    else
        tnext = t_prev;
    
end