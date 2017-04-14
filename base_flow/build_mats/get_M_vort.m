function M_vort = get_M_vort( grid_parms )

%get a diagonal matrix M_vort that converts vorticity to circulation via
%                           omega = M_vort * gamma
%(i.e. M_vort has 1/h^2 terms on the diagonal)

%Inputs: grid_parms -- data structure containing m (number of points in x
%dirn), n (number of points in y dirn), mg (number of grid levels), and len
%(length of domain in x-dirn ==> dx = len / m )

m = grid_parms.m; n = grid_parms.n; mg = grid_parms.mg; len = grid_parms.len;

del = len / m;

%Get size of M_vort
nrows = get_vort_ind( m, n-1, mg, grid_parms );
ncols = nrows ;

M_vort = sparse( nrows, ncols );

for glev = 1 : mg
    
    %grid spacing on current grid
    delb = del * 2.d0^(glev-1);
    
    if glev == 1
        ind_s = 1;
    else
        ind_s = 1 + get_vort_ind(m-1,n-1,glev-1,grid_parms);
    end
    
    if glev < mg
        ind_e = get_vort_ind(m-1,n-1,glev,grid_parms);
    else
        ind_e = get_vort_ind(m,n-1,glev,grid_parms);
    end
    
    ind = ind_s : ind_e;
    
    M_vort = M_vort + sparse( ind, ind, (1/delb)^2 * ones(size(ind)), nrows, ncols);
    
end





