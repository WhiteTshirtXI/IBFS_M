function C = get_C_main( C,  grid_parms )

m = grid_parms.m; n = grid_parms.n; mg = grid_parms.mg;
nrows = grid_parms.nrows;
ncols = grid_parms.ncols;


%--First build block corresponding to x-velocities

    %vorticity points above x-velocity point
    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    velx_ind = get_velx_ind(indvelx,indvely, 1, grid_parms);

    indvortx = repmat(1 : m-1,[1,n-1]);
    indvorty = repelem(1 : n-1, m-1);
    vort_ind = get_vort_ind(indvortx,indvorty, 1, grid_parms);

    C = C + sparse( velx_ind, vort_ind, ...
        ones(size(vort_ind)), nrows, ncols);

    %vorticity points below x-velocity point
    indvelx = repmat(1 : m-1, [1, n-1]);
    indvely = repelem(2 : n, m-1);
    velx_ind = get_velx_ind(indvelx,indvely, 1 , grid_parms);

    indvortx = repmat(1 : m-1, [1, n-1]);
    indvorty = repelem(1 : n-1, m-1);
    vort_ind = get_vort_ind(indvortx,indvorty,1,grid_parms);

    C = C - sparse( velx_ind, vort_ind,...
        ones(size(vort_ind)), nrows, ncols);

%--    

%--Now build y-velocity block

    %rows start at end of x-velocity block:
    n_add = get_velx_ind( m-1, n, 1, grid_parms );

    %vorticity points to the right of y-velocity point
    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    vely_ind = get_vely_ind(indvelx,indvely,1,grid_parms);

    indvortx = repmat(1 : m-1,[1,n-1]);
    indvorty = repelem(1 : n-1, m-1);
    vort_ind = get_vort_ind(indvortx,indvorty,1,grid_parms);

    C = C - sparse(n_add + vely_ind, vort_ind, ...
        ones(size(vort_ind)), nrows, ncols);

    %vorticity points to the left of y-velocity point
    indvelx = repmat(2 : m,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    vely_ind = get_vely_ind(indvelx,indvely,1,grid_parms);

    indvortx = repmat(1 : m-1,[1,n-1]);
    indvorty = repelem(1 : n-1, m-1);
    vort_ind = get_vort_ind(indvortx,indvorty,1,grid_parms);

    C = C + sparse(n_add + vely_ind, vort_ind, ...
        ones(size(vort_ind)), nrows, ncols);

%--

