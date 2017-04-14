function Q = get_Q_fineedge_BCs( Q, glev, grid_parms )

%BC's for Q from coarser grid needed on current grid


m = grid_parms.m; n = grid_parms.n; mg = grid_parms.mg;
nrows = grid_parms.nrows;
ncols = grid_parms.ncols;

%y-velocity index starts at the end of all x-velocities
n_add = get_velx_ind( m-1, n, mg, grid_parms );

%--x-velocity block needs top and bottom BCs
    %Bottom part
    
        %get coarse grid y-velocity indices on bottom:
        indvelx = m/4+1 : 3*m/4;
        indvely = n/4;
        vely_ind = get_vely_ind(indvelx,indvely,glev+1,grid_parms);
    
        %points that need to two points from coarser domain
        indvelx = 2 : 2 : m-2;
        indvely = 1;
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms); 

        Q = Q - 1/4 * sparse( velx_ind, n_add + vely_ind(1:end-1), ...
                ones(size(velx_ind)), nrows, ncols ) ...
              - 1/4 * sparse( velx_ind, n_add + vely_ind(2:end), ...
                ones(size(velx_ind)), nrows, ncols ); 
        
        %Points that just need one point from coarser
        indvelx = 1 : 2 : m-1;
        indvely = 1;
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms); 

        Q = Q - 1/2 * sparse( velx_ind, n_add + vely_ind, ...
        ones(size(velx_ind)), nrows, ncols );

    %Top part
    
        %get coarse grid y-velocity indices on bottom:
        indvelx = m/4+1 : 3*m/4;
        indvely = 3*n/4;
        vely_ind = get_vely_ind(indvelx,indvely,glev+1,grid_parms);
    
        %points that need to two points from coarser domain
        indvelx = 2 : 2 : m-2;
        indvely = n;
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms); 

        Q = Q - 1/4 * sparse( velx_ind, n_add + vely_ind(1:end-1), ...
                ones(size(velx_ind)), nrows, ncols ) ...
              - 1/4 * sparse( velx_ind, n_add + vely_ind(2:end), ...
                ones(size(velx_ind)), nrows, ncols ); 
        
        %Points that just need one point from coarser
        indvelx = 1 : 2 : m-1;
        indvely = n;
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms); 

        Q = Q - 1/2 * sparse( velx_ind, n_add + vely_ind, ...
        ones(size(velx_ind)), nrows, ncols );

%--

%--y-velocity block needs left and right BCs

    %Left part
    
        %get coarse grid x-velocity indices on left:
        indvelx = m/4;
        indvely = n/4+1 : 3*n/4;
        velx_ind = get_velx_ind(indvelx,indvely,glev+1,grid_parms);
    
        %points that need to two points from coarser domain
        indvelx = 1;
        indvely = 2 : 2 : n-2;
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms); 

        Q = Q + 1/4 * sparse( n_add + vely_ind, velx_ind(1:end-1), ...
                ones(size(vely_ind)), nrows, ncols ) ...
              + 1/4 * sparse( n_add + vely_ind, velx_ind(2:end), ...
                ones(size(vely_ind)), nrows, ncols ); 
        
        %Points that just need one point from coarser
        indvelx = 1;
        indvely = 1 : 2 : n-1;
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms); 

        Q = Q + 1/2 * sparse( n_add + vely_ind, velx_ind, ...
        ones(size(velx_ind)), nrows, ncols );
    
    %Right part
    
        %get coarse grid x-velocity indices on right:
        indvelx = 3*m/4;
        indvely = n/4+1 : 3*n/4;
        velx_ind = get_velx_ind(indvelx,indvely,glev+1,grid_parms);
    
        %points that need to two points from coarser domain
        indvelx = m;
        indvely = 2 : 2 : n-2;
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms); 

        Q = Q + 1/4 * sparse( n_add + vely_ind, velx_ind(1:end-1), ...
                ones(size(vely_ind)), nrows, ncols ) ...
              + 1/4 * sparse( n_add + vely_ind, velx_ind(2:end), ...
                ones(size(vely_ind)), nrows, ncols ); 
        
        %Points that just need one point from coarser
        indvelx = m;
        indvely = 1 : 2 : n-1;
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms); 

        Q = Q + 1/2 * sparse( n_add + vely_ind, velx_ind, ...
        ones(size(velx_ind)), nrows, ncols );


%--

