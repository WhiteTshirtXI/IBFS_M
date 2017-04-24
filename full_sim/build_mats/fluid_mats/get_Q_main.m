function Q = get_Q_main( Q, grid_parms )

%Build part of Q (without accounting for BC's). 
%Note that grid_parms contains info on m, n, mg, and number of
%rows and columns of Q.

m = grid_parms.m; n = grid_parms.n; 
nrows = grid_parms.nrows;
ncols = grid_parms.ncols;

%y-velocity index starts at the end of all x-velocities
n_add = get_velx_ind( m-1, n, 1, grid_parms );


%--Block corresponding to x-velocities

    %y-vels to bottom left of current x-vel point
    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(2 : n, m-1);
    velx_ind = get_velx_ind(indvelx,indvely,1,grid_parms);

    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    vely_ind = get_vely_ind(indvelx,indvely,1,grid_parms);

    Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
        ones(size(vely_ind)), nrows, ncols);

    %y-vels to bottom right of current x-vel point
    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(2 : n, m-1);
    velx_ind = get_velx_ind(indvelx,indvely,1,grid_parms);

    indvelx = repmat(2 : m,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    vely_ind = get_vely_ind(indvelx,indvely,1,grid_parms);

    Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
        ones(size(vely_ind)), nrows, ncols);

    %y-vels to top left of current x-vel point
    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    velx_ind = get_velx_ind(indvelx,indvely,1,grid_parms);

    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    vely_ind = get_vely_ind(indvelx,indvely,1,grid_parms);

    Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
        ones(size(vely_ind)), nrows, ncols);

    %y-vels to top right of current x-vel point
    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    velx_ind = get_velx_ind(indvelx,indvely,1,grid_parms);

    indvelx = repmat(2 : m,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    vely_ind = get_vely_ind(indvelx,indvely,1,grid_parms);

    Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
        ones(size(vely_ind)), nrows, ncols);

%--

%--Block corresponding to y-velocities

    %x-vels to bottom left of current y-vel point
    indvelx = repmat(2 : m,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    vely_ind = get_vely_ind(indvelx,indvely,1,grid_parms);

    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    velx_ind = get_velx_ind(indvelx,indvely,1,grid_parms);

    Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
        ones(size(vely_ind)), nrows, ncols);

    %x-vels to bottom right of current y-vel point
    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    vely_ind = get_vely_ind(indvelx,indvely,1,grid_parms);

    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    velx_ind = get_velx_ind(indvelx,indvely,1,grid_parms);

    Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
        ones(size(vely_ind)), nrows, ncols);

    %x-vels to top left of current y-vel point
    indvelx = repmat(2 : m,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    vely_ind = get_vely_ind(indvelx,indvely,1,grid_parms);

    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(2 : n, m-1);
    velx_ind = get_velx_ind(indvelx,indvely,1,grid_parms);

    Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
        ones(size(vely_ind)), nrows, ncols);

    %x-vels to top right of current y-vel point
    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(1 : n-1, m-1);
    vely_ind = get_vely_ind(indvelx,indvely,1,grid_parms);

    indvelx = repmat(1 : m-1,[1,n-1]);
    indvely = repelem(2 : n, m-1);
    velx_ind = get_velx_ind(indvelx,indvely,1,grid_parms);

    Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
        ones(size(vely_ind)), nrows, ncols);

%--

