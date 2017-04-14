function velx_ind = get_velx_ind( i, j, glev, grid_parms )

m = grid_parms.m;
n = grid_parms.n;

%Return vector indices for x - velocity (velx_ind) 
%corresponding to 2d array indices. 

% e.g. get_velx_ind(2, 3, 1, grid_parms) returns 2*(m-1) + 2, which is the 
% vector index of the x velocity at the point with xindex = 2, yindex = 3,
% and gridlevel = 1.

%Note: x velocity on grid 1 is of size (m-1)*n.
%   For all other grids it is of size (m-1)*n/2 + m/2*n/2.


%1st gridlevel is easy...
if glev == 1
    
      velx_ind = (m-1) .* (j-1) + i ;
    
%If not 1st grid level...
else
    
    %index starts from the end of the last gridlevel...
    n_add = (m-1) * n + ... %contribution from 1st gridlevel
        (glev - 2) * ( (m-1)*n/2 + m/2*n/2 ); %contribution from remaining glevs 
    
    
    %points that are below finer grid:
    velx_ind = (j <= n/4) .* ( n_add + (m-1) .* (j-1) + i );    
    
    %points that contain finer grid and are to left of finer grid:
    n_bott = (m-1)*n/4;
    velx_ind = velx_ind + ( j > n/4 & j <= 3*n/4 & i <= m/4 ) .* ...
        ( n_add + n_bott + m/2*(j-n/4-1) + i );
    
    %points that contain finer grid and are to right of finer grid:
    velx_ind = velx_ind + (j > n/4 & j <= 3*n/4 & i >= 3*m/4) .* ...
        ( n_add + n_bott + m/2*(j-n/4-1) + (i-m/2+1) );
    
    %points that are above finer grid:
    n_mid = m/2 * n/2;
    velx_ind = velx_ind + (j > 3*n/4 ).* ...
        ( n_add + n_bott + n_mid + (j-3*n/4-1)*(m-1) + i );
    
    
end


