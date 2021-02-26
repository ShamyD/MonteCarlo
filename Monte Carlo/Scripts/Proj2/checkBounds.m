function within_bounds = checkBounds(array, lb, rb)
    logical_lower = array>=lb;
    logical_higher = array<=rb;
    within_bounds = logical_lower+logical_higher;
    within_bounds = within_bounds==2;
    
end

