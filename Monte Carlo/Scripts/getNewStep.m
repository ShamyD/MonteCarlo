function step = getNewStep(N,d)

    s = [1 0;0 1;-1 0;0 -1];

    rI = randi(4, N, 1);
    step = s(rI, :);
end

