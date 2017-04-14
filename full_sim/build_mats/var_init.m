function [beta0, L0] = var_init( x, y )

%Update geometry based on original coordinates (x,y) in global frame and
%generalized displacement vector u:

nel = length(x) - 1;

beta0 = zeros( nel, 1);
L0 = zeros( nel, 1);

for j = 1 : nel
    
    dx = x(j+1) - x(j);
    dy = y(j+1) - y(j);
    
    L0(j) = sqrt( dx^2 + dy^2 );
    
    beta0( j ) = atan2( dy, dx );
    
end