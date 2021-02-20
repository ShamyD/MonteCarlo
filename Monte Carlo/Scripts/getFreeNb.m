function [free_nb, nr_free_nb] = getFreeNb(walks, N, d, dirs)
    
    
    curr_p = walks(:,end-d+1:end);
    prev_ps = walks(:,1:end-d);
    
    %Get nb:s
    nbs = repmat(curr_p,1,d*2) + repmat(dirs,N,1);
    size_dirs = size(dirs);
    size_prev_ps = size(prev_ps);
    
    % Matrix with with 1 if dir is okay and 0 otherwise.
    poss_dir = zeros(N,d*d*2);
    
    
    %Loop over neighbours.
    for dir=1:(size_dirs(2))/d
        
        nb_2_check = nbs(:,d*dir - d + 1:d*dir);
        nb_is_same = zeros(N,(size_prev_ps(2))/d);
        for col=1:(size_prev_ps(2))/d
            old_col = prev_ps(:, col*d -d + 1:col*d);
            temp_col = old_col - nb_2_check;
            members = ismember(temp_col,zeros(1,d), 'rows');
            
            nb_is_same(:,col) = members;
        end
        
        nb_ok_temp = sum(nb_is_same',1)';
        
        %nb_ok has 1s when the point is unvisited and zeros otherwise.
        nb_ok = nb_ok_temp==0;
        
        poss_dir(:,dir*d-d+1:dir*d) = repmat(nb_ok,1,d);        
    end
    nr_free_nb = sum(poss_dir')'./d;
    
    poss_dir(poss_dir == 0) = NaN;
    %Matrix to store free neighbours in.
    free_nb = poss_dir.*nbs;
    
    
end


