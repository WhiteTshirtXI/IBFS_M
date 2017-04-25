function grid_parms = setup_reg( xb, grid_parms )

%Return the regularization matrix ET (scaled to be the transpose of E),
%which takes quantities on the IB and smears them to the flow domain.


%parameters used throughout function
m = grid_parms.m; n = grid_parms.n; nb = grid_parms.nb; nf = 2 * nb;
len = grid_parms.len; supp = grid_parms.supp; offx = grid_parms.offx;
offy = grid_parms.offy; del = len / m; supp = ceil(supp / 2); 

%index of flow grid points nearest to IB points
indexx = zeros( nf, 1 );
indexx( 1 : nb ) = round( (xb(1:nb) + offx) ./ del );
indexx( nb + 1 : nf ) = round( ( xb( nb + 1 : nf) + offy ) ./del ) ;

grid_parms.indexx = indexx;

%build weights
grid_parms.wght = zeros( (2*supp + 1)*(2*supp + 1), nf );
for i = 1 : nb
    
    nextx = 0;
    nexty = 0;
    
    for j = -supp : supp

        for k = -supp : supp

            %--rows corresponding to x-vels
                x = del * ( indexx(i)-1 + k ) - offx; %x posn of grid point
                y = del * ( indexx(i + nb)-1 + j ) - offy - del/2; %y posn 

                nextx = nextx + 1;
                %entries to put into columns
                grid_parms.wght( nextx, i)= delta_h( x, xb(i), del) .* ...
                    delta_h( y, xb(i + nb), del) ;
            %--
            
            %--rows corresponding to y-vels
                x = del * ( indexx(i)-1 + k ) - offx - del/2; %x posn 
                y = del * ( indexx(i + nb)-1 + j ) - offy; %y posn 

                %entries to put into columns
                nexty = nexty + 1;
                grid_parms.wght( nexty, i + nb ) = ...
                    delta_h( x, xb(i), del) .* ...
                    delta_h( y, xb(i + nb), del) ;  
            %--
        end

    end

end

