function regTx = get_regT( x, grid_parms )

%Return E * x, where E is the interpolation matrix that takes quantities
%from the flow domain and interpolates it onto the IB.

%parameters used throughout function
m = grid_parms.m; n = grid_parms.n; nb = grid_parms.nb; nf = 2 * nb;
supp = grid_parms.supp; supp = ceil(supp / 2); 

regTx = zeros( nf, 1);
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
            regTx( k ) = regTx( k ) + wght( next, k ) * x( indvelx );
            
            indvely = get_vely_ind( i+p-1, j+el-1, 1, grid_parms ) + ...
                get_velx_ind( m-1, n, 1, grid_parms );
            regTx( k+nb ) = regTx( k+nb ) + wght( next, k+nb )...
                *x( indvely );
            
        end

    end

end

