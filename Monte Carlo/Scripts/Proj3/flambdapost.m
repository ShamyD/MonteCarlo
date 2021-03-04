function lambda_post = flambdapost(d, theta, t, ni)
    % Sample from the posterior distribution of lambda given
    % observations of the rest.
    % ni is a function of tau.
    
    lambda_post = zeros(1,d);
    
    for i =1:d
        lambda_post(i) = gamrnd(ni(i) + 2,1/(t(i+1)-t(i) + theta));
    end
end

