function [wave_height] = gumb_inv(mu, beta, s)
    % Evaluates the inverse of the gumbel cdf.
    
    wave_height = mu - beta.*log(-log(s));
    
end

