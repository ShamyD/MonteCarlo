function c_cum = SISR_sample(N,n,d)

% Generate N walks from g(x_(0:n))
walks = zeros(N,d*(n + 1));
%Store weights.
weight = ones(N,n+1);

dirs = getDirs(d);

    for c = 1:n
        
        % Get which free nb:s there are and how many (for all walks).
        [free_nb, nr_free_nb] = getFreeNb(walks(:,1:d*c),N,d,dirs);
    
        % Get new point for each of the N separate walks.
        prev_col = walks(:,d*c - d + 1:d*c);
        new_col = getFreeStep(N,d,prev_col, free_nb, nr_free_nb);
        walks(:,d*c + 1:d*c + d) = new_col;
        
        % calculate weights
        weight(:,c+1) = nr_free_nb;
    
        %Resample walks with help of weights.
        mult = resampling(weight(:, c));
        matr = multToMatr(mult);
        walks = matr*walks;
    
    end
    
    % Average weights for each length between 1:n.
    result = sum(weight,1)/N;

    % Create c_sisr sample for n between 1:n.
    c_cum = ones(1, n + 1);
    c_cum(1,1) = result(1,1);
    for i = 2:n + 1
        c_cum(1,i) = c_cum(1,i-1)*result(1,i);
    end

    c_cum = c_cum(2:end);
end

