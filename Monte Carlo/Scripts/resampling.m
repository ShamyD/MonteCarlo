function multiplicity = resampling(weights)
%   Function that takes in weights and uses them to resample
%   based on said weights.
    
    N = length(weights);
    sum_weights = sum(weights,1);
    %Vector with discrete distribution.
    prob = weights/sum_weights;
    
    %Generate prob_cum which is the cumulative distribution with
    % 1 as it's last element. Done using triangular matrix.
    tri = tril(ones(N,N));
    prob_cum = tri*prob;
    %Repeat said matrix many times.
    prob_cum_rep = repmat(prob_cum, 1, N);
    
    %Generate samples from uniform distribution.
    randvec = rand(1,N);
    % Repeat many times.
    randmatr = repmat(randvec, N,1);
    
    % Subtract. The number of new negative elements on each row is
    % the number of copies of that walk that survives.
    diff_matr = randmatr - prob_cum_rep;
    cum_samp = sum(diff_matr < 0, 2);
    
    %Subtract the number of negative elements on previous row from
    % nr of negative elements on current row to get nr of copies
    %that survive resampling.
    distributor = diag(ones(N,1)) + diag(-1*ones(N-1,1),-1);
    multiplicity = distributor*cum_samp;
   
end

