function [free_steps] = getFreeStep(N,d,prev_points, free_nb, nr_free_nb)
    
    
    which_nb = ceil(rand(N,1).*nr_free_nb);
    which_nb = (which_nb*2 - 1);
    
    % free_steps = free_nb_t([which_nb which_nb+1]);
    
    
    for row=1:N
        r = free_nb(row,:);
        r = [r(~isnan(r)) r(isnan(r))];
        
        nr_points = nr_free_nb(row);
        if nr_points ~= 0
            free_steps(row,:) = r(which_nb(row):which_nb(row) + 1);
        else
            free_steps(row,:) = prev_points(row, :);
        end
    end 
end

