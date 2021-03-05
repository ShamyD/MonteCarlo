function [t_posterior] = ftpost(t, lambda, ni)
    % Evaluate the logarithm of the posterior of t given lambda,
    % and theta.
    
    endpoints_interval = t(2:end);
    starts_interval = t(1:end-1);
    
    interval_sizes = endpoints_interval - starts_interval;
    lambda_interval_sizes = lambda .* interval_sizes;
    sum_lambda_interval_sizes = sum(lambda_interval_sizes);
    

    prod_lambdai_to_ni = sum(ni.*log(lambda));
    
    if issorted(t)
        ft_prior = sum(log((interval_sizes)));
    else
        ft_prior = log(0);
    end
    
    t_posterior = -sum_lambda_interval_sizes + prod_lambdai_to_ni + ft_prior;
end

