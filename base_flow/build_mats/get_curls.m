function [C, R] = get_curls( grid_parms )

%Build C, R
%C -- discrete curl matrix that takes curl of the streamfunction and maps 
%     it onto velocity edges.
%R -- discrete curl matrix that takes curl of velocity flux to get
%     circulation
%     (R = C' except for BCs that need to be accounted for).

%Inputs: grid_parms -- data structure containing m (number of points in x
%dirn), n (number of points in y dirn), mg (number of grid levels), and len
%(length of domain in x-dirn ==> dx = len / m )

m = grid_parms.m; n = grid_parms.n; mg = grid_parms.mg;

%Get size of C
nrows = get_velx_ind( m-1, n, mg, grid_parms ) + ...
    get_vely_ind( m, n-1, mg, grid_parms ) ;
ncols = get_vort_ind( m-1, n-1, mg, grid_parms ) ;

grid_parms.nrows = nrows;
grid_parms.ncols = ncols;

C = sparse( nrows, ncols );

%C matrix without accounting for BC's
%--Note: can't incorporate the BC's in this for loop because the transpose
%        Cp uses different BC's...
for glev = 1 : mg
    
    %Main blocks from current grid (without BC's)
    C = get_C_main( C, glev, grid_parms );
    
end

%Get the transpose operator (also a curl operator) 
R = C';

%Modify C and R to account for BCs
for glev = 1 : mg
   
    %BC's for C from coarser grid needed on current grid
    if mg > 1 & glev < mg
        
        C = get_C_fineedge_BCs( C, glev, grid_parms );

    end
    
    %BC's for R from finer grid needed on current grid
    if glev > 1
        
        R = get_R_coarseinner_BCs( R, glev, grid_parms );

    end
    
    
end



grid_parms = rmfield( grid_parms, 'nrows'); 
grid_parms = rmfield( grid_parms, 'ncols');
