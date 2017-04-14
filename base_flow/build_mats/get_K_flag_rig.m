function K = get_K_flag_rig( R_E, R_th, R_sh, L0, c, s, L, qL , parms)

%Build global stiffness matrix from current beam configuration...

%This builds the DIMENSIONLESS stiffness matrix. The dimensional form has
%EA and EI as factors that it gets multiplied by. In the dimensionless form
%of the equations, EA gets replaced by R_E*R_th ( E/(rho_f*U^2) ) and EI gets
%replaced by R_E*R_sh ( I/(L^3*b) )

nel = length( c );
K = zeros( 3*(nel + 1) );
for j = 1 : nel
    
    r_c = R_sh/R_th;
    CL = R_E*R_th/L0 * [2 0 0; 0 4*r_c 2*r_c; 0 2*r_c 4*r_c];
    B = [-c(j) -s(j) 0 c(j) s(j) 0; ...
        -s(j)/L(j) c(j)/L(j) 1 s(j)/L(j) -c(j)/L(j) 0;...
        -s(j)/L(j) c(j)/L(j) 0 s(j)/L(j) -c(j)/L(j) 1];
    
    K1 = B' * CL * B;
    
    qind = 3*( j - 1 ) + 1 : 3*(j-1) + 3;
    Nf = qL( qind( 1 ) );
    Mf1 = qL( qind( 2 ) );
    Mf2 = qL( qind( 3 ) );
    
    zvect = [s(j) -c(j) 0 -s(j) c(j) 0]';
    rvect = -[c(j)  s(j) 0 -c(j) -s(j) 0]';
    
    K2 = Nf/L(j)*(zvect*zvect') + (Mf1 + Mf2)/( L(j)^2 )*( rvect*zvect' + ...
        zvect*rvect' );
    
    K_e = K1 + K2;
    
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
            K(:,j) = 0;
            K(j,j) = 1;
        end
        
    end
    
else
    
    for j = 1 : 3
       
        count = count + 1;
        
        if bc_type(count) == 1
            K(j,:) = 0;
            K(:,j) = 0;
            K(j,j) = 1;
        end
        
    end
    
    
end


