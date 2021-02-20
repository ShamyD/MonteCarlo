function [cum_sums] = SISR_sampling(N, n, d, rep)
    cum_sums = zeros(rep, n);

    for i = 1:rep
        cum_sums(i, :) = SISR_sample(N, n, d);
        if mod(i, 10) == 0
            i
        end
    end
end

