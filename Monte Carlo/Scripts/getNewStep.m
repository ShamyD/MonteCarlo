function step = getNewStep(N,d)
    % Gets a new step uniformly from all 4 neighbours for all walks.
    s = [1 0;0 1;-1 0;0 -1];

    rI = randi(4, N, 1);
    step = s(rI, :);
end

