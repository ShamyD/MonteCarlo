function theta = fthetapost(lambda, phi)
    % Sample from theta posterior, given observed values of lamba, t, tau.
    d = size(lambda,2);
    theta = gamrnd(2+2*d,phi + sum(lambda));    
end

