function [free_nb, nr_free_nb] = getFreeNb(walks, N, d, dirs)
    %Determines which neighbours are free for each of the latest
    %points in all walks. 
    %walks - matrix of the walk up to this point. (N x d*c) in size
    %dirs - matrix with all possible directions. Depends on 
    %dimension.
    
    
    %curr_p is coordinates for latest point
    curr_p = walks(:,end-d+1:end);
    % Previous points.
    prev_ps = walks(:,1:end-d);
    
    %Get nb:s for current points.
    nbs = repmat(curr_p,1,d*2) + repmat(dirs,N,1);
    
    size_dirs = size(dirs);
    size_prev_ps = size(prev_ps);
    
    % Matrix with with 1 if dir is okay and 0 otherwise. It's
    % d 1:s if okay and d zeros if not.
    poss_dir = zeros(N,d*d*2);
    
    
    %Loop over all directions and determine for each neighbour
    %direction if free. 
    for dir=1:(size_dirs(2))/d
        
        
        nb_2_check = nbs(:,d*dir - d + 1:d*dir);
        nb_is_same = zeros(N,(size_prev_ps(2))/d);
        %Loop over previous points and determine if any 
        %is the same as nb_2_check. If so nb_is_same gets
        %a 1 otherwise 0.
        for col=1:(size_prev_ps(2))/d
            old_col = prev_ps(:, col*d -d + 1:col*d);
            temp_col = old_col - nb_2_check;
            members = ismember(temp_col,zeros(1,d), 'rows');
            
            nb_is_same(:,col) = members;
        end
        
        %Sum over rows, and if not zero the neighbour is not free.
        nb_ok_temp = sum(nb_is_same',1)';
        
        %nb_ok has 1s when the point is unvisited and zeros otherwise.
        nb_ok = nb_ok_temp==0;
        
        poss_dir(:,dir*d-d+1:dir*d) = repmat(nb_ok,1,d);        
    end
    % Sum for number of free neighbours.
    nr_free_nb = sum(poss_dir')'./d;
    
    poss_dir(poss_dir == 0) = NaN;
    %Matrix to store free neighbours in.
    free_nb = poss_dir.*nbs;
    
    
end


