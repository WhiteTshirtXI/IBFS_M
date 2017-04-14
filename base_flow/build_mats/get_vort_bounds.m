function vort_bounds = get_vort_bounds( glev, grid_parms )

%Get indices on coarse grid level (glev) that corresponds to the left,
%right, top, and bottom indices on the finer grid level.

m = grid_parms.m; n = grid_parms.n;

vort_bounds.left = get_vort_ind( m/4, n/4 : 3*n/4, glev, grid_parms );
vort_bounds.right = get_vort_ind( 3*m/4, n/4 : 3*n/4, glev, grid_parms );
vort_bounds.bottom = get_vort_ind( m/4 : 3*m/4, n/4, glev, grid_parms );
vort_bounds.top = get_vort_ind( m/4 : 3*m/4, 3*n/4, glev, grid_parms );





