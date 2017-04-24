function [C, R] = get_curls( grid_parms )

%Build C, R
%C -- discrete curl matrix that takes curl of the streamfunction and maps 
%     it onto velocity edges.
%R -- discrete curl matrix that takes curl of velocity flux to get
%     circulation
%     (R = C').

%Inputs: grid_parms -- data structure containing m (number of points in x
%dirn), n (number of points in y dirn), mg (number of grid levels), and len
%(length of domain in x-dirn ==> dx = len / m )

m = grid_parms.m; n = grid_parms.n; 

%Get size of C
nrows = get_velx_ind( m-1, n, 1, grid_parms ) + ...
    get_vely_ind( m, n-1, 1, grid_parms ) ;
ncols = get_vort_ind( m-1, n-1, 1, grid_parms ) ;

grid_parms.nrows = nrows;
grid_parms.ncols = ncols;

C = sparse( nrows, ncols );

%Main blocks from current grid (without BC's)
C = get_C_main( C, grid_parms );
 
%Get the transpose operator (also a curl operator) 
R = C';

