function M = get_M_flag( R_rho, R_th, h0, nel, nb, parms )


%Builds the mass matrix for the DIMENSIONLESS system of equations.
%The dimensional form has a factor rho*A that pre-multiplies the matrix,
%where rho is the mass per unit volume and A is the cross sectional area.

%In the dimensionless form of the equations, this factor is replaced with 
%R_rho*R_th, where R_rho is the ratio of the solid density to the fluid
%density (rho_s/rho_f) and R_th is the thickness ratio h/L (h is the
%thickness of the beam and L is the length of the beam)

M = zeros(3*(nel + 1), 3*(nel + 1) );
for j = 1 : nel
    
    M_e = R_rho*R_th*h0/420 * ...
       [140 0 0 70 0 0; ...
        0 156 22*h0 0 54 -13*h0; ...
        0  22*h0 4*h0^2 0 13*h0 -3*h0^2; ...
        70  0  0 140 0 0; ...
        0  54  13*h0  0 156 -22*h0; ...
        0  -13*h0  -3*h0^2  0  -22*h0  4*h0^2];
    
    index = get_ind_flag( j );
    
    M = M_flag_assemble( M, M_e, index );
    
end

%Apply BCs:


%clamped?
if parms.clamped == 'T'
   
    bc_type = [1 1 1];
    
else
    
    bc_type = [1 1 0];
    
end

count = 0;
if parms.inverted == 'T'
   
    for j = 3*nb - 2 : 3*nb
       
        count = count + 1;
        
        if bc_type(count) == 1
            M(j,:) = 0;
            M(:,j) = 0;
        end
        
    end
    
else
    
    for j = 1 : 3
       
        count = count + 1;
        
        if bc_type(count) == 1
            M(j,:) = 0;
            M(:,j) = 0;
        end
        
    end
    
    
end


