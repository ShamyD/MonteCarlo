function [free_steps] = getFreeStep(N,d,prev_points, free_nb, nr_free_nb)
    
    
    which_nb = ceil(rand(N,1).*nr_free_nb);
    which_nb = (which_nb*d - d + 1);
    
    free_steps = zeros(N,d);
    
    
    for row=1:N
        r = free_nb(row,:);
        r = [r(~isnan(r)) r(isnan(r))];
        
        nr_points = nr_free_nb(row);
        if nr_points ~= 0
            free_steps(row,:) = r(which_nb(row):which_nb(row) + d - 1);
        else
            free_steps(row,:) = prev_points(row, :);
        end
    end 
end

