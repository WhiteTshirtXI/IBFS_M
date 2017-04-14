function C = get_C_fineedge_BCs( C, glev, grid_parms )

%BC's for C from coarser grid needed on current grid


m = grid_parms.m; n = grid_parms.n; mg = grid_parms.mg;
nrows = grid_parms.nrows;
ncols = grid_parms.ncols;

%Define vorticity indices on coarse grid that correspond to bottom, left, 
%right, and top edges of finer grid.
vort_bounds = get_vort_bounds( glev + 1, grid_parms );

%--x-velocity block
    %Bottom part
        %points that need to average coarser domain
        indvelx = 1 : 2 : m-1;
        indvely = 1;
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);                

        C = C - 1/2*sparse(velx_ind, vort_bounds.bottom(1:end-1), ...
            ones(size(velx_ind)), nrows, ncols) - ...
            1/2*sparse(velx_ind, vort_bounds.bottom(2:end), ...
            ones(size(velx_ind)), nrows, ncols);

        %Points that don't need to be averaged:
        indvelx = 2 : 2 : m-2;
        indvely = 1;
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

        C = C - sparse(velx_ind, vort_bounds.bottom(2:end-1), ...
            ones(size(velx_ind)), nrows, ncols);


    %Top part
        %points that need to average coarser domain
        indvelx = 1 : 2 : m-1;
        indvely = n;
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

        C = C + 1/2*sparse(velx_ind, vort_bounds.top(1:end-1), ...
            ones(size(velx_ind)), nrows, ncols) + ...
            1/2*sparse(velx_ind, vort_bounds.top(2:end), ...
            ones(size(velx_ind)), nrows, ncols);

        %Points that don't need to be averaged:
        indvelx = 2 : 2 : m-2;
        indvely = n;
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

        C = C + sparse(velx_ind, vort_bounds.top(2:end-1), ...
            ones(size(velx_ind)), nrows, ncols);
        
 %--
 
 %--y-velocity block
 
    %rows start at end of x-velocity block:
        n_add = get_velx_ind( m-1, n, mg, grid_parms );
 
    %Left part
        %points that need to average coarser domain
        indvelx = 1;
        indvely = 1 : 2 : n-1;
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);                

        C = C + 1/2*sparse(n_add + vely_ind, vort_bounds.left(1:end-1), ...
            ones(size(vely_ind)), nrows, ncols) + ...
            1/2*sparse(n_add + vely_ind, vort_bounds.left(2:end), ...
            ones(size(vely_ind)), nrows, ncols);

        %Points that don't need to be averaged:
        indvelx = 1;
        indvely = 2 : 2 : n-2;
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

        C = C + sparse(n_add + vely_ind, vort_bounds.left(2:end-1), ...
            ones(size(vely_ind)), nrows, ncols);
 
        
    %Right part
        %points that need to average coarser domain
        indvelx = m;
        indvely = 1 : 2 : n-1;
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);                

        C = C - 1/2*sparse(n_add + vely_ind, vort_bounds.right(1:end-1), ...
            ones(size(vely_ind)), nrows, ncols) - ...
            1/2*sparse(n_add + vely_ind, vort_bounds.right(2:end), ...
            ones(size(vely_ind)), nrows, ncols);

        %Points that don't need to be averaged:
        indvelx = m;
        indvely = 2 : 2 : n-2;
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

        C = C - sparse(n_add + vely_ind, vort_bounds.right(2:end-1), ...
            ones(size(vely_ind)), nrows, ncols);
 
 

 %--
 
 
        
        