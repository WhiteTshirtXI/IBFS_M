function M = M_assemble( M, M_e, index )

for i = 1 : length(index)
   
    ii = index(i);
    
    for j = 1 : length(index)
       
        jj = index(j);
        
        M( ii, jj ) = M( ii, jj ) + M_e( i, j );
        
    end
    
end
