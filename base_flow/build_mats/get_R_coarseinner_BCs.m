function R = get_R_coarseinner_BCs( R, glev, grid_parms )

%Get BC's from finer grid on interior edges of current grid.


m = grid_parms.m; n = grid_parms.n; mg = grid_parms.mg;
nrows = grid_parms.nrows;
ncols = grid_parms.ncols;
n_add = get_velx_ind( m-1, n, mg, grid_parms );


%Define vorticity indices on current grid that correspond to bottom, left, 
%right, and top edges of finer grid.
vort_bounds = get_vort_bounds( glev, grid_parms );

%--
%Left edge (vorticity on edge requires y velocity from finer domain)
indvelx = repmat(1:2,[1, n/2-1]);
indvely = repelem(2 : 2 : n-2, 2);
vely_ind = get_vely_ind(indvelx,indvely,glev-1,grid_parms);

vort_ind = repelem(vort_bounds.left(2:end-1), 2);

%note: nrows and ncols are switched because we are working with R, not C.
R = R + 1/2 * sparse( vort_ind, n_add + vely_ind, ...
    ones(size(vely_ind)), ncols, nrows );


indvelx = repmat(1:2,[1, n/2-1]);
indvely = repelem(1 : 2 : n-3, 2);
vely_ind = get_vely_ind(indvelx,indvely,glev-1,grid_parms);

vort_ind = repelem(vort_bounds.left(2:end-1), 2);

%note: nrows and ncols are switched because we are working with R, not C.
R = R + 1/4 * sparse( vort_ind, n_add + vely_ind, ...
    ones(size(vely_ind)), ncols, nrows );


indvelx = repmat(1:2,[1, n/2-1]);
indvely = repelem(3 : 2 : n-1, 2);
vely_ind = get_vely_ind(indvelx,indvely,glev-1,grid_parms);

vort_ind = repelem(vort_bounds.left(2:end-1), 2);

%note: nrows and ncols are switched because we are working with R, not C.
R = R + 1/4 * sparse( vort_ind, n_add + vely_ind, ...
    ones(size(vely_ind)), ncols, nrows );

%--

%--
%Right edge (vorticity on edge requires y velocity from finer domain)
indvelx = repmat(m-1: m,[1, n/2-1]);
indvely = repelem(2 : 2 : n-2, 2);
vely_ind = get_vely_ind(indvelx,indvely,glev-1,grid_parms);

vort_ind = repelem(vort_bounds.right(2:end-1), 2);

%note: nrows and ncols are switched because we are working with R, not C.
R = R - 1/2 * sparse( vort_ind, n_add + vely_ind, ...
    ones(size(vely_ind)), ncols, nrows );


indvelx = repmat(m-1: m,[1, n/2-1]);
indvely = repelem(1 : 2 : n-3, 2);
vely_ind = get_vely_ind(indvelx,indvely,glev-1,grid_parms);

vort_ind = repelem(vort_bounds.right(2:end-1), 2);

%note: nrows and ncols are switched because we are working with R, not C.
R = R - 1/4 * sparse( vort_ind, n_add + vely_ind, ...
    ones(size(vely_ind)), ncols, nrows );


indvelx = repmat(m-1: m,[1, n/2-1]);
indvely = repelem(3 : 2 : n-1, 2);
vely_ind = get_vely_ind(indvelx,indvely,glev-1,grid_parms);

vort_ind = repelem(vort_bounds.right(2:end-1), 2);

%note: nrows and ncols are switched because we are working with R, not C.
R = R - 1/4 * sparse( vort_ind, n_add + vely_ind, ...
    ones(size(vely_ind)), ncols, nrows );


%--

%--
%Bottom edge (vorticity on edge requires x velocity from finer domain)
indvelx = repmat(2 : 2 : m-2,[1, 2]);
indvely = repelem(1 : 2, m/2-1);
velx_ind = get_velx_ind(indvelx,indvely,glev-1,grid_parms);

vort_ind = repmat(vort_bounds.bottom(2:end-1), [1,2]);

%note: nrows and ncols are switched because we are working with R, not C.
R = R - 1/2 * sparse( vort_ind, velx_ind, ...
    ones(size(velx_ind)), ncols, nrows );

indvelx = repmat(1 : 2 : m-3,[1, 2]);
indvely = repelem(1 : 2, m/2-1);
velx_ind = get_velx_ind(indvelx,indvely,glev-1,grid_parms);

vort_ind = repmat(vort_bounds.bottom(2:end-1), [1,2]);

%note: nrows and ncols are switched because we are working with R, not C.
R = R - 1/4 * sparse( vort_ind, velx_ind, ...
    ones(size(velx_ind)), ncols, nrows );

indvelx = repmat(3 : 2 : m-1,[1, 2]);
indvely = repelem(1 : 2, m/2-1);
velx_ind = get_velx_ind(indvelx,indvely,glev-1,grid_parms);

vort_ind = repmat(vort_bounds.bottom(2:end-1), [1,2]);

%note: nrows and ncols are switched because we are working with R, not C.
R = R - 1/4 * sparse( vort_ind, velx_ind, ...
    ones(size(velx_ind)), ncols, nrows );

%--

%--
%Top edge (vorticity on edge requires x velocity from finer domain)
indvelx = repmat(2 : 2 : m-2,[1, 2]);
indvely = repelem(n-1 : n, m/2-1);
velx_ind = get_velx_ind(indvelx,indvely,glev-1,grid_parms);

vort_ind = repmat(vort_bounds.top(2:end-1), [1,2]);

%note: nrows and ncols are switched because we are working with R, not C.
R = R + 1/2 * sparse( vort_ind, velx_ind, ...
    ones(size(velx_ind)), ncols, nrows );


indvelx = repmat(1 : 2 : m-3,[1, 2]);
indvely = repelem(n-1 : n, m/2-1);
velx_ind = get_velx_ind(indvelx,indvely,glev-1,grid_parms);

vort_ind = repmat(vort_bounds.top(2:end-1), [1,2]);

%note: nrows and ncols are switched because we are working with R, not C.
R = R + 1/4 * sparse( vort_ind, velx_ind, ...
    ones(size(velx_ind)), ncols, nrows );


indvelx = repmat(3 : 2 : m-1,[1, 2]);
indvely = repelem(n-1 : n, m/2-1);
velx_ind = get_velx_ind(indvelx,indvely,glev-1,grid_parms);

vort_ind = repmat(vort_bounds.top(2:end-1), [1,2]);

%note: nrows and ncols are switched because we are working with R, not C.
R = R + 1/4 * sparse( vort_ind, velx_ind, ...
    ones(size(velx_ind)), ncols, nrows );
%--





