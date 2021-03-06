function W = get_W_main( W, grid_parms )

%Build part of W (without accounting for BC's). 
%Note that grid_parms contains info on m, n, mg, and number of
%rows and columns of W.

m = grid_parms.m; n = grid_parms.n; 
nrows = grid_parms.nrows;
ncols = grid_parms.ncols;

%y-velocity index starts at the end of all x-velocities
n_add = get_velx_ind( m-1, n, 1, grid_parms );


%-- x-vel block

    %vorticity points above x-velocity point
    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    velx_ind = get_velx_ind(indvelx,indvely,1,grid_parms);

    indvortx = repmat(1 : m-1,[1,n-1]);
    indvorty = repelem(1 : n-1, m-1);
    vort_ind = get_vort_ind(indvortx,indvorty,1,grid_parms);

    W = W + 1/2*sparse( velx_ind, vort_ind, ...
        ones(size(vort_ind)), nrows, ncols);

    %vorticity points below x-velocity point
    indvelx = repmat(1 : m-1, [1, n-1]);
    indvely = repelem(2 : n, m-1);
    velx_ind = get_velx_ind(indvelx,indvely,1,grid_parms);

    indvortx = repmat(1 : m-1, [1, n-1]);
    indvorty = repelem(1 : n-1, m-1);
    vort_ind = get_vort_ind(indvortx,indvorty,1,grid_parms);

    W = W + 1/2*sparse( velx_ind, vort_ind,...
        ones(size(vort_ind)), nrows, ncols);

%--

%-- y-vel block

    %vorticity points to the right of y-velocity point
    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    vely_ind = get_vely_ind(indvelx,indvely,1,grid_parms);

    indvortx = repmat(1 : m-1,[1,n-1]);
    indvorty = repelem(1 : n-1, m-1);
    vort_ind = get_vort_ind(indvortx,indvorty,1,grid_parms);

    W = W + 1/2*sparse(n_add + vely_ind, vort_ind, ...
        ones(size(vort_ind)), nrows, ncols);

    %vorticity points to the left of y-velocity point
    indvelx = repmat(2 : m,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    vely_ind = get_vely_ind(indvelx,indvely,1,grid_parms);

    indvortx = repmat(1 : m-1,[1,n-1]);
    indvorty = repelem(1 : n-1, m-1);
    vort_ind = get_vort_ind(indvortx,indvorty,1,grid_parms);

    W = W + 1/2*sparse(n_add + vely_ind, vort_ind, ...
        ones(size(vort_ind)), nrows, ncols);

%--

