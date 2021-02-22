function [free_steps] = getFreeStep(N,d,prev_points, free_nb, nr_free_nb)
    %Takes arguments:
    %prev_points - The last points from which to take next step.
    %free_nb - Matrix containing the coordinates of each free
    %neighbour on rows. If not free - NaN.
    
    %Randomize index for which neighbour to be chosen.
    which_nb = ceil(rand(N,1).*nr_free_nb);
    which_nb = (which_nb*d - d + 1);
    
    %Storage matrix.
    free_steps = zeros(N,d);
    
    %Loop over walks.
    for row=1:N
        %Make row r have non NaN values first. 
        r = free_nb(row,:);
        r = [r(~isnan(r)) r(isnan(r))];
        
        nr_points = nr_free_nb(row);
        %If nr of free nb:s is zero set new point to last point.
        %Else choose point with random generated indexes.
        if nr_points ~= 0
            free_steps(row,:) = r(which_nb(row):which_nb(row) + d - 1);
        else
            free_steps(row,:) = prev_points(row, :);
        end
    end 
end

