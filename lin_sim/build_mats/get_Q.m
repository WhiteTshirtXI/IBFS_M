function Q = get_Q( grid_parms )

%Build averaging operator Q = [Qx; Qy]. The x-velocity block (Qx) takes y
%velocities and averages them onto the x-velocity edges, and the y-velocity
%block (Qy) takes x velocities and averages them onto y-velocity edges.

%Note: this does not give the Q used in the code. The Q used in the code is
%postmultiplies Minv to the Q obtained here.

%Inputs: grid_parms -- data structure containing m (number of points in x
%dirn), n (number of points in y dirn), and mg (number of grid levels)

m = grid_parms.m; n = grid_parms.n; mg = grid_parms.mg;

%Get size of Q
nrows = get_velx_ind( m-1, n, mg, grid_parms ) + ...
    get_vely_ind( m, n-1, mg, grid_parms );
ncols = nrows;

grid_parms.nrows = nrows;
grid_parms.ncols = ncols;

%Initialize
Q = sparse( nrows, ncols );

%Loop through gridlevels
for glev = 1 : mg
    
    %Get main blocks (not accounting for BC's)
    Q = get_Q_main( Q, glev, grid_parms );
    
    %BC contributions from coarser grid on current grid
    if mg > 1 & glev < mg
        
        Q = get_Q_fineedge_BCs( Q, glev, grid_parms );
        
    end
    
    %BC contributions from finer grid on current grid
    if glev > 1
        
        Q = get_Q_coarseinner_BCs( Q, glev, grid_parms );

    end
end


grid_parms = rmfield( grid_parms, 'nrows'); 
grid_parms = rmfield( grid_parms, 'ncols');