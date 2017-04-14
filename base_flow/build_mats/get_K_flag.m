function K = get_K_flag( E, A, I, c, s, L, parms )



%Build global stiffness matrix from current beam configuration...

nel = length( c );
K = zeros( 3*(nel + 1) );

for j = 1 : nel
    
    K1 = zeros(6,6);
    K1(1,1) = E*A/L(j);
    K1(1,4) = - E*A/L(j);
    
    K1(2,2) = 12*E*I/(L(j)^3);
    K1(2,3) = 6*E*I/(L(j)^2);
    K1(2,5) = -12*E*I/(L(j)^3);
    K1(2,6) = 6*E*I/(L(j)^2);
    
    K1(3,2) = K1(2,3);
    K1(3,3) = 4*E*I/L(j);
    K1(3,5) = -6*E*I/(L(j)^2);
    K1(3,6) = 2*E*I/L(j);
    
    K1(4,1) = K1(1,4);
    K1(4,4) = E*A/L(j);
    
    K1(5,2) = K1(2,5);
    K1(5,3) = K1(3,5);
    K1(5,5) = 12*E*I/(L(j)^3);
    K1(5,6) = -6*E*I/(L(j)^2);
    
    K1(6,2) = K1(2,6);
    K1(6,3) = K1(3,6);
    K1(6,5) = K1(5,6);
    K1(6,6) = 4*E*I/L(j);


    R = [c(j) s(j) 0 0 0 0; ...
         -s(j) c(j) 0 0 0 0; ...
         0 0 1 0 0 0; ...
         0 0 0 c(j) s(j) 0; ...
         0 0 0 -s(j) c(j) 0; ...
         0 0 0 0 0 1];

    
    K_e = R'* K1 * R;
    
    %Assemble the stiffness matrix:
    index = get_ind_flag( j );
    K = K_flag_assemble( K, K_e, index );
    
    
    
end

%BCs

%clamped?
if parms.clamped == 'T'
   
    bc_type = [1 1 1];
    
else
    
    bc_type = [1 1 0];
    
end

nb = nel+1;
count = 0;

if parms.inverted == 'T'
   
    for j = 3*nb - 2 : 3*nb
       
        count = count + 1;
        
        if bc_type(count) == 1
            K(j,:) = 0;
%             K(:,j) = 0;
            K(j,j) = 1;
        end
        
    end
    
else
    
    for j = 1 : 3
       
        count = count + 1;
        
        if bc_type(count) == 1
            K(j,:) = 0;
%             K(:,j) = 0;
            K(j,j) = 1;
        end
        
    end
    
    
end


