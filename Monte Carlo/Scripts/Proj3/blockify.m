function [blockedMatrix] = blockify(matrix,nbrBlocks)
    
    s = size(matrix);
    blockedMatrix = zeros(nbrBlocks, s(2));
    lenB = floor(s(1)/nbrBlocks);
    
    for i = 1:nbrBlocks
        top = i*lenB;
        bot = (i-1)*lenB+1;
        if i == nbrBlocks
            top = s(1);
        end
        m = mean(matrix(bot:top,:),1);
        blockedMatrix(i, :) = m;    
    end
    
    
end

