function ET = get_ET2( xb, grid_parms )

%Return the regularization matrix ET (scaled to be the transpose of E),
%which takes quantities on the IB and smears them to the flow domain.


%parameters used throughout function
m = grid_parms.m; n = grid_parms.n; nb = grid_parms.nb; nf = 2 * nb;
len = grid_parms.len; supp = grid_parms.supp; offx = grid_parms.offx;
offy = grid_parms.offy; del = len / m; supp = ceil(supp / 2); 


%Get size of ET
nrows = get_velx_ind( m-1, n, 1, grid_parms ) + ...
    get_vely_ind( m, n-1, 1, grid_parms );
ncols = length(xb);

ET = sparse( nrows, ncols );

%index of flow grid points nearest to IB points
indexx = zeros( nf, 1 );
indexx( 1 : nb ) = round( (xb(1:nb) + offx) ./ del );
indexx( nb + 1 : nf ) = round( ( xb( nb + 1 : nf) + offy ) ./del ) ;


%build ET
for i = 1 : nb
    
    for j = -supp : supp

        for k = -supp : supp

            %--rows corresponding to x-vels
            
                x = del * ( indexx(i)-1 + k ) - offx; %x posn of grid point
                y = del * ( indexx(i + nb)-1 + j ) - offy - del/2; %y posn 

                
                %entries to put into columns
                del_h = delta_h( x, xb(i), del) .* ...
                    delta_h( y, xb(i + nb), del) ;
                
                %Add to ET:
                indvel = get_velx_ind( indexx(i)-1 + k, ...
                    indexx(i + nb)-1 + j, 1, grid_parms );
                ET = ET + sparse( indvel, i, del_h, nrows, ncols );
            %--
            
            %--rows corresponding to y-vels
                x = del * ( indexx(i)-1 + k ) - offx - del/2; %x posn 
                y = del * ( indexx(i + nb)-1 + j ) - offy; %y posn 

                %entries to put into columns
                del_h = delta_h( x, xb(i), del) .* ...
                    delta_h( y, xb(i + nb), del) ;

                %Add to ET:
                %y index starts after all the x-vels
                n_add = get_velx_ind( m-1, n, 1, grid_parms );
                
                indvel = n_add + get_vely_ind( indexx(i)-1 + k, ...
                    indexx(i + nb)-1 + j, 1, grid_parms );
                ET = ET + sparse( indvel, i + nb, del_h, nrows, ncols );
            %--
        end

    end

end

