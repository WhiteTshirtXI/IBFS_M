function vort_ind = get_vortb_ind( i, j, glev, grid_parms )

m = grid_parms.m;
n = grid_parms.n;
mg = grid_parms.mg;

%Return vector indices for vorticity (vort_ind) 
%corresponding to 2d array indices. 

% e.g. get_vort_ind(2, 3, 1, grid_parms) returns 2*(m-1) + 2, which is the 
% vector index of the x velocity at the point with xindex = 2, yindex = 3,
% and gridlevel = 1.

%Note: vorticity on grid 1 is of size (m-1)*(n-1).
%   For all other grids it is of size (m-1)*n/2 + m/2*(n/2-1).



vort_ind = m .* (j-1) + i ;
