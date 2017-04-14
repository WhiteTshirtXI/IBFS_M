function M_vel = get_M_vel( grid_parms )

%get a diagonal matrix M_vel that converts velocity to flux via
%                           u = M * q
% (i.e. M has 1/h terms)

%Inputs: grid_parms -- data structure containing m (number of points in x
%dirn), n (number of points in y dirn), mg (number of grid levels), and len
%(length of domain in x-dirn ==> dx = len / m )

m = grid_parms.m; n = grid_parms.n; mg = grid_parms.mg; len = grid_parms.len;

del = len / m;

%Get size of C
nrows = get_velx_ind( m-1, n, mg, grid_parms ) + ...
    get_vely_ind( m, n-1, mg, grid_parms ) ;
ncols = nrows ;

M_vel = sparse( nrows, ncols );

for glev = 1 : mg
    
    %grid spacing on current grid
    delb = del * 2.d0^(glev-1);
    
    %-- x-vel block
        if glev == 1
            ind_s = 1;
        else
            ind_s = 1 + get_velx_ind(m-1,n,glev-1,grid_parms);
        end

        ind_e = get_velx_ind(m-1,n,glev,grid_parms);
        
        ind = ind_s : ind_e;
        
        M_vel = M_vel + sparse( ind, ind, (1/delb) * ones(size(ind)), nrows, ncols);
    %--
    
    %-- y-vel block
        if glev == 1
            ind_s = 1;
        else
            ind_s = 1 + get_vely_ind(m,n-1,glev-1,grid_parms);
        end

        ind_e = get_vely_ind(m,n-1,glev,grid_parms);
        
        ind = ind_s : ind_e;
        
        %y-vel index starts after x-vels
        n_add = get_velx_ind(m-1,n,mg, grid_parms);
        ind = ind + n_add;
        
        M_vel = M_vel + sparse( ind, ind, (1/delb) * ones(size(ind)), nrows, ncols);
    %--
    
    
    
end





