function dirs = getDirs(d)

    dirs=zeros(1,2*d^2);

    for i=0:d-1
        dirs(2*i*d+1:2*i*d+d) = [zeros(1,i) 1 zeros(1,d-i-1)];
        dirs(2*i*d+d+1:2*i*d+d+d) = -1*[zeros(1,i) 1 zeros(1,d-i-1)];
    end    
end

