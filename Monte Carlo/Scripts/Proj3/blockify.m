function [blockedMatrix] = blockify(matrix,nbrBlocks)
    
    s = size(matrix);
    blockedMatrix = zeros(nbrBlocks, s(2));
    iter = ceil(s(1)/nbrBlocks)
    
    for i = 1:iter
        top = i*nbrBlocks;
        m = mean(matrix(,:),1)
    end
    
    
end

