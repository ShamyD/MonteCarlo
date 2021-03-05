function ni = create_ni(tau,t)
    % Create vector with nr of accidents in each interval.
    nr_of_accidents = size(tau,2);
    d = size(t,2) - 1;
    ends_of_intervals = t(2:end);
    nr_of_intervals = size(ends_of_intervals,2);
    
    t_prim = ends_of_intervals';
    rep_t = repmat(t_prim, 1, nr_of_accidents);
    rep_tau = repmat(tau,nr_of_intervals,1);
    
    % non_zero_matr is zero if year of accident is smaller than 
    non_zero_matr = rep_tau-rep_t;
    non_zero_matr = non_zero_matr < 0;
    accs_cumulative = sum(non_zero_matr,2);
    
    % accs_cumulative now has nr of accidents in that interval and
    % those earlier than that. Multiply by matrix to make
    % the cumulative nr of accidents into accidents per interval.
    double_diag = diag(ones(1,d)) - diag(ones(1,d-1),-1);
    ni = (double_diag * accs_cumulative)';
    
end

