function multiplicity = resampling(weights)
%RESAMPLING Summary of this function goes here
%   Detailed explanation goes here
    N = length(weights);
    sum_weights = sum(weights,1);
    prob = weights/sum_weights;
    
    tri = tril(ones(N,N));
    prob_cum = tri*prob;
    prob_cum_rep = repmat(prob_cum, 1, N);
    
    randvec = rand(1,N);
    randmatr = repmat(randvec, N,1);
    
    diff_matr = randmatr - prob_cum_rep;
    cum_samp = sum(diff_matr < 0, 2);
    
    distributor = diag(ones(N,1)) + diag(-1*ones(N-1,1),-1);
    multiplicity = distributor*cum_samp;
   
end

