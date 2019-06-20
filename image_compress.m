function [ ] = image_compress( image_filename, destination_filename )
%image_compress compress the file "image_filename" into "destination_filename"

    % Opening file
    B = double(imread(image_filename));
    B = B(50:249,50:249); % Comment this line if you want to compress the whole image
    nr = size(B,1);
    nc = size(B,2);
    
    %%%%%%%%%%%%%%%%%
    % Decorrelating %
    E = zeros(nr,nc);
    Evertcat = [];
    for j=1:nc
        for i=1:nr
            causalNeighborhood = [];
            if i > 1
                causalNeighborhood = [causalNeighborhood, B(i-1, j)];
            end
            if j > 1
                causalNeighborhood = [causalNeighborhood, B(i, j-1)];
            end
            if i > 1 && j > 1
                causalNeighborhood = [causalNeighborhood, B(i-1, j) + B(i, j-1) - B(i-1, j-1)];
            end
            if i == 1 && j == 1
                causalNeighborhood = [0];
            end
            E(i,j) = B(i,j) - median(causalNeighborhood);
        end
        % Concatening every column
        Evertcat = vertcat(Evertcat, E(:,j));
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Looking for best L, with best p
    L = (50:50:2000);
    CL = [];
    
    best_size = 1e10;
    best_L = 1;
    
    for block_size_idx = 1:size(L, 2)
        block_size = L(block_size_idx);
        nb_blocks = ceil(size(Evertcat,1) / block_size);
        total_size = 0;
        optimals_p = [];
        
        for i = 1:nb_blocks
            % Last block must be troncated
            block = Evertcat(1 + (i-1)*block_size : min(i*block_size, 200*200));
            optimal_block_size = 1e10;
            optimal_block_p = 1;
            % Looking for optimal p %
            for p = 1:16
                bsize = GRestimation(block, p);
                if bsize < optimal_block_size
                    optimal_block_size = bsize;
                    optimal_block_p = p;
                end
            end
            optimals_p = [optimals_p, optimal_block_p];
            total_size = total_size + optimal_block_size;
        end
        % Final size is sum of optimal block sizes + encoding optimal p %
        % Each p value can be transmitted using the number of bits needed to encode the maximum value p %
        p_encoding_size = floor(log2(max(optimals_p))) + 1;
        final_size = total_size + nb_blocks * p_encoding_size;
        CL = [CL, final_size];
        if final_size < best_size 
            best_size = final_size;
            best_p = optimals_p;
            best_L = block_size;
        end
    end
    
    %%%%%%%%%%%%%%%
    % Compression %
    
    compressed_bits = [];
    nb_blocks = ceil(size(Evertcat,1) / best_L);
    
    for i = 1:nb_blocks
        p = best_p(i);
        block = Evertcat(1 + (i-1)*best_L : min(i*best_L, 200*200));
        
        for intgr_idx = 1:length(block)
            intgr = block(intgr_idx);
            if intgr < 0; sign = 1; intgr = -intgr; else; sign = 0; end
            least_significant = mod(intgr, 2^p);
            least_significant_binary = fliplr(de2bi(least_significant));
            most_significant = floor(intgr / 2^p);
            most_significant_unary = [zeros(1, most_significant), 1];
            
            compressed_bits = [compressed_bits, sign];
            n = 0;
            while n < (p - length(least_significant_binary))
                compressed_bits = [compressed_bits, 0];
                n = n + 1;
            end
            compressed_bits = horzcat(compressed_bits, least_significant_binary);
            compressed_bits = horzcat(compressed_bits, most_significant_unary);
        end
    end
    
    %%%%%%%%%%%%%%%
    % Saving file %
    
    fid = fopen(destination_filename, 'w');
    fwrite(fid, compressed_bits, 'ubit1');
    fclose(fid);
end