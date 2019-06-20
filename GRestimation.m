function [ rsize ] = GRestimation( block, p )
%GRestimation Compute the size the block of integrer would take with parameter p
    rsize = 0;
    for i = 1:size(block,1)
        intgr = block(i);
        rsize = rsize + 1 + p + (floor(abs(intgr)/(2^p)) + 1);
    end
end