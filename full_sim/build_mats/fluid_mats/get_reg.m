function regx = get_reg( x, grid_parms )

%Return ET * x, where ET is the regularization matrix that takes quantities
%from the IB and smears them onto the flow domain.

%parameters used throughout function
m = grid_parms.m; n = grid_parms.n; nb = grid_parms.nb; nf = 2 * nb;
len = grid_parms.len; supp = grid_parms.supp; supp = ceil(supp / 2); 
nq = (m-1)*n + (n-1) * m;

regx = zeros( nq, 1);
indexx = grid_parms.indexx; 
wght = grid_parms.wght;

for k = 1 : nb
    
    i = indexx( k );
    j = indexx( k + nb );
    
    next = 0;
    
    for el = -supp : supp

        for p = -supp : supp

            next = next + 1;
            indvelx = get_velx_ind( i+p-1, j+el-1, 1, grid_parms );
            regx( indvelx ) = regx( indvelx ) + wght( next, k ) * x( k );
            
            indvely = get_vely_ind( i+p-1, j+el-1, 1, grid_parms ) + ...
                get_velx_ind( m-1, n, 1, grid_parms );
            regx( indvely ) = regx( indvely ) + wght( next, k+nb )*x( k+nb );
            
        end

    end

end

