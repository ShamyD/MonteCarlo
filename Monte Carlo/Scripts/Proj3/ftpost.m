function [t_posterior] = ftpost(t, lambda, ni)
    '-------'
    % Evaluate posterior of t given lambda, theta and theta.
    
    endpoints_interval = t(2:end);
    starts_interval = t(1:end-1);
    
    interval_sizes = endpoints_interval - starts_interval;
    lambda_interval_sizes = lambda .* interval_sizes;
    sum_lambda_interval_sizes = sum(lambda_interval_sizes);
    
    prod_lambdai_to_ni = prod(lambda .^ni);
    
    ft_prior = prod(interval_sizes);
    
    t_posterior = exp(-sum_lambda_interval_sizes)*prod_lambdai_to_ni*ft_prior;
end

