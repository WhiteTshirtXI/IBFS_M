function K = K_flag_assemble( K, K_e, index )

for i = 1 : length(index)
   
    ii = index(i);
    
    for j = 1 : length(index)
       
        jj = index(j);
        
        K( ii, jj ) = K( ii, jj ) + K_e( i, j );
        
    end
    
end
