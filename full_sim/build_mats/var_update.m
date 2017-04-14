function [c, s, L, qL, Fint] = var_update(x, y, u, R_E, R_th, R_sh, ...
        beta0, L0, parms )

%Update geometry based on original coordinates (x,y) in global frame and
%generalized displacement vector u:

%This subroutine works in dimensionless parameters. The dimensional form of
%this subroutine uses the terms E*A and E*I. In the dimensionless version
%E*A is replaced with R_E*R_th and E*I is replaced with R_E*R_sh (see get_M 
%subroutine for a description of these dimensionless parameters)

nel = length(x) - 1;
nb = parms.nb;
c = zeros( nel, 1);
s = zeros( nel, 1);
L = zeros( nel, 1);
qL = zeros( 3*nel, 1);
Fint = zeros(3*(nel + 1), 1); 

for j = 1 : nel
    
    index = get_ind_flag(j);
    
    dx = x(j+1) + u(index(4)) - ( x(j) + u(index(1)) );
    dy = y(j+1) + u(index(5)) - ( y(j) + u(index(2)) );
    
    L(j) = sqrt( dx^2 + dy^2 );
    
    c(j) = dx/L(j);
    s(j) = dy/L(j);
    
    uL = ( L(j)^2 - L0^2 )/( L(j) + L0 );
    
    qind = 3*(j-1)+1 : 3*(j-1)+3;
    qL(qind(1)) = R_E*R_th*uL/( L0 );
    
    %             beta = atan2( dy, dx );
    beta1 = u(index(3)) + beta0( j );
    beta2 = u(index(6)) + beta0( j );
    theta1L = atan2( c(j)*sin(beta1) - s(j)*cos(beta1), ...
        c(j)*cos(beta1) + s(j)*sin(beta1) );
    theta2L = atan2( c(j)*sin(beta2) - s(j)*cos(beta2), ...
        c(j)*cos(beta2) + s(j)*sin(beta2) );
    
    qL(qind(2)) = 2*R_E*R_sh/L0 * (2*theta1L + theta2L);
    qL(qind(3)) = 2*R_E*R_sh/L0 * (theta1L + 2*theta2L);
    
    B = [-c(j) -s(j) 0 c(j) s(j) 0; ...
        -s(j)/L(j) c(j)/L(j) 1 s(j)/L(j) -c(j)/L(j) 0;...
        -s(j)/L(j) c(j)/L(j) 0 s(j)/L(j) -c(j)/L(j) 1];
    
    qint = B'*qL(qind(1):qind(3) ); %Internal forces in global frame
            
    Fint = F_assemble( Fint, qint, index );
    
end

%apply BCs to Fint

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
            Fint(j) = 0;
        end

    end

else

    for j = 1 : 3

        count = count + 1;

        if bc_type(count) == 1
            Fint(j) = 0;
        end

    end


end


