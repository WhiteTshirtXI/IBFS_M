function Q = get_Q_coarseinner_BCs( Q, glev, grid_parms )

%BC's for Q from finer grid needed on current grid


m = grid_parms.m; n = grid_parms.n; mg = grid_parms.mg;
nrows = grid_parms.nrows;
ncols = grid_parms.ncols;

%y-velocity index starts at the end of all x-velocities
n_add = get_velx_ind( m-1, n, mg, grid_parms );


%-- x-velocity block needs left and right BCs

    %left edge
    
        %get fine grid y-vel indices on left:
        indvelx = repmat(1 : 2, [1, n/2-1]);
        indvely = repelem(2 : 2 : n-2, 2);
        vely_ind = n_add + ...
            get_vely_ind( indvelx, indvely, glev-1, grid_parms );
        
        %get current grid x-vel indices on left
        indvelx = m/4;
        indvely = n/4+1 : 3*n/4;
        velx_ind = get_velx_ind( indvelx, indvely, glev, grid_parms );
        
        %y-vel contribution on bottom right of x-vel point
        velx_use = repelem( velx_ind(2:end), 2);
        
        Q = Q - 1/8 * sparse( velx_use, vely_ind, ...
            ones(size(velx_use)), nrows, ncols );
        
        %y-vel contribution on top right of x-vel point
        velx_use = repelem( velx_ind(1:end-1), 2);
        
        Q = Q - 1/8 * sparse( velx_use, vely_ind, ...
            ones(size(velx_use)), nrows, ncols );
        
        
    %right edge
    
        %get fine grid y-vel indices on right:
        indvelx = repmat(m-1 : m, [1, n/2-1]);
        indvely = repelem(2 : 2 : n-2, 2);
        vely_ind = n_add + ...
            get_vely_ind( indvelx, indvely, glev-1, grid_parms );
        
        %get current grid x-vel indices on right
        indvelx = 3*m/4;
        indvely = n/4+1 : 3*n/4;
        velx_ind = get_velx_ind( indvelx, indvely, glev, grid_parms );
        
        %y-vel contribution on bottom left of x-vel point
        velx_use = repelem( velx_ind(2:end), 2);
        
        Q = Q - 1/8 * sparse( velx_use, vely_ind, ...
            ones(size(velx_use)), nrows, ncols );
        
        %y-vel contribution on top left of x-vel point
        velx_use = repelem( velx_ind(1:end-1), 2);
        
        Q = Q - 1/8 * sparse( velx_use, vely_ind, ...
            ones(size(velx_use)), nrows, ncols );
        
%--

%-- y-velocity block needs top and bottom BCs

    %bottom edge
    
        %get fine grid x-vel indices on bottom:
        indvelx = repmat(2 : 2 : m-2, [1, 2]);
        indvely = repelem(1 : 2, m/2-1);
        velx_ind = get_velx_ind( indvelx, indvely, glev-1, grid_parms );
        
        %get current grid y-vel indices on bottom
        indvelx = m/4 + 1 : 3*m/4;
        indvely = n/4;
        vely_ind = n_add + ...
            get_vely_ind( indvelx, indvely, glev, grid_parms );
        
        %x-vel contribution on top left of y-vel point
        vely_use = repmat( vely_ind(2:end), [1,2]);
        
        Q = Q + 1/8 * sparse( vely_use, velx_ind, ...
            ones(size(vely_use)), nrows, ncols );
        
        %x-vel contribution on top right of y-vel point
        vely_use = repmat( vely_ind(1:end-1), [1,2]);
        
        Q = Q + 1/8 * sparse( vely_use, velx_ind, ...
            ones(size(vely_use)), nrows, ncols );
        
        
    %top edge
    
        %get fine grid x-vel indices on top:
        indvelx = repmat(2 : 2 : m-2, [1, 2]);
        indvely = repelem(n-1 : n, m/2-1);
        velx_ind = get_velx_ind( indvelx, indvely, glev-1, grid_parms );
        
        %get current grid y-vel indices on top
        indvelx = m/4 +1 : 3*m/4;
        indvely = 3*n/4;
        vely_ind = n_add + ...
            get_vely_ind( indvelx, indvely, glev, grid_parms );
        
        %x-vel contribution on bottom left of y-vel point
        vely_use = repmat( vely_ind(2:end), [1,2]);
        
        Q = Q + 1/8 * sparse( vely_use, velx_ind, ...
            ones(size(vely_use)), nrows, ncols );
        
        %x-vel contribution on bottom right of y-vel point
        vely_use = repmat( vely_ind(1:end-1), [1,2]);
        
        Q = Q + 1/8 * sparse( vely_use, velx_ind, ...
            ones(size(vely_use)), nrows, ncols );


%--



