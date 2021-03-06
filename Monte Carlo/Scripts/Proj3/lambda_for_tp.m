function lambda_t = lambda_for_tp(lambda, d, t, tp)
    % Takes a lambda vector, breakpoint vector t, and a timepoint
    % vector tp and evaluates the lambda in each point.
%     size(lambda)
%     size(d)
%     size(t)
%     size(tp)
    
    
    inner_tp = tp(2:end-1);
    l_tp = size(inner_tp,2);
    lambda_index = zeros(1,l_tp);
    
    for i=2:d+1
        
        lambda_index = lambda_index + double(inner_tp<t(i));
    end
    lambda_index = d - lambda_index + 1;
    lambda_t = lambda(lambda_index);
end

