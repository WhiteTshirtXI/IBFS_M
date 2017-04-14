function F = F_assemble( F, F_e, index )

for i = 1 : length(index)
   
    ii = index(i);
      
    F( ii ) = F( ii ) + F_e( i );
     
end        