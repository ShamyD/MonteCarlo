function c_cum = SISR_sample(N,n)

% Generate N walks from g(x_(0:n))
d = 2;
walks = zeros(N,d*(n + 1));
weight = ones(N,n+1);
result = ones(1,n);

    for c = 1:n
        
        % Get which free nb:s there are and how many (for all walks).
        [free_nb, nr_free_nb] = getFreeNb(walks(:,1:d*c),N,d);
    
        % Get new point for each of the N separate walks.
        prev_col = walks(:,d*c - d + 1:d*c);
        new_col = getFreeStep(N,d,prev_col, free_nb, nr_free_nb);
        walks(:,d*c + 1:d*c + d) = new_col;
        
        % calculate weights
        %N_tot(:, c) = nr_free_nb;
%       if c>1
        weight(:,c+1) = nr_free_nb;
%         end
    
    
%         if c ~= 0
            mult = resampling(weight(:, c));
            matr = multToMatr(mult);
            walks = matr*walks;
%         end
    
   
 
    
    end

    result = sum(weight,1)/N;

    c_cum = ones(1, n + 1);
    c_cum(1,1) = result(1,1);
    for i = 2:n + 1
        c_cum(1,i) = c_cum(1,i-1)*result(1,i);
    end

    c_cum = c_cum(2:end);
end

