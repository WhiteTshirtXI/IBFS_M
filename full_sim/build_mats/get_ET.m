function ET = get_ET( grid_parms )

%Return the regularization matrix ET (scaled to be the transpose of E),
%which takes quantities on the IB and smears them to the flow domain.
m = grid_parms.m; n = grid_parms.n; 
nb = grid_parms.nb;
len = grid_parms.len;
xb = grid_parms.xb;
supp = grid_parms.supp; offx = grid_parms.offx;
offy = grid_parms.offy;

del = len / m;

%Get size of ET
nrows = get_velx_ind( m-1, n, 1, grid_parms ) + ...
    get_vely_ind( m, n-1, 1, grid_parms );
ncols = length(xb);

ET = sparse( nrows, ncols );

%--rows corresponding to x-vels

    %x and y points on physical grid (for fluid domain)
    xu = (del : del : (m-1)*del ) - offx;
    yu = (del/2 : del : (n-1/2) * del ) - offy;
    
    %x and y points on physical grid (for IB)
    xb_x = xb( 1 : nb );
    xb_y = xb( 1 + nb : 2*nb );
    
    %find x points within support
    xu_r = repmat( xu, [1,length(xb_x)]);
    xb_rx = repelem( xb_x, length(xu) );
    ind_supp_x = ( abs( xu_r - xb_rx ) <= supp * del );
    
    %Get index of these points for fluid grid
    i_supp = round( (xu_r(ind_supp_x) + offx) / del );
    
    %Get index of these points for IB:
    i_supp_xbx = round( ind_supp_x .* repelem( 1:nb, length(xu) ) );
    i_supp_xbx = i_supp_xbx( i_supp_xbx ~= 0 );
        
    %find y points within support
    yu_r = repmat( yu, [1,length(xb_y)]);
    xb_ry = repelem( xb_y, length(yu) );
    ind_supp_y = ( abs( yu_r - xb_ry ) <= supp * del );
    
    %Get y-index of these points for fluid grid
    j_supp = round( (yu_r(ind_supp_y) + offy) /del + 1/2  );
                    
    %Get index of these points for IB:
    j_supp_xby = round( ind_supp_y .* repelem( 1:nb, length(yu) ) );
    j_supp_xby = j_supp_xby( j_supp_xby ~= 0 );
    
    %For each IB point, add nonzero weights...
    for j = 1 : nb
        
        %x-indices on IB corresponding to current body point
        ind_xbx = (i_supp_xbx == j);
        
        %x-indices on flow grid that are within support of IB point
        ind_x = i_supp( ind_xbx );
        
        %y-indices on IB corresponding to current body point
        ind_xby = (j_supp_xby == j);
        
        %y-indices on flow grid that are within support of IB point
        ind_y = j_supp( ind_xby );
        
        %Combine flow indices
        indvelx = repmat( ind_x, [1,length(ind_y)]);
        indvely = repelem( ind_y, length(ind_x) );
        
        %rows of ET to add to
        velx_ind = get_velx_ind( indvelx, indvely, 1, grid_parms );
        
        %columns to add to
        xb_ind = j * ones( size(velx_ind ) );
        
        %entries to put into columns
        del_h = delta_h( xu(indvelx), xb_x(xb_ind), del) .* ...
            delta_h( yu(indvely), xb_y(xb_ind), del) ;
        
        
        %Add to ET:
        %y index starts after all the x-vels
        ET = ET + sparse( velx_ind, xb_ind, del_h, nrows, ncols );

        
    end
%--

%--rows corresponding to y-vels

    %x and y points on physical grid (for fluid domain)
    xv = (del/2 : del : (m-1/2) * del ) - offx;
    yv = (del : del : (n-1)*del ) - offy;
    
    %x and y points on physical grid (for IB)
    xb_x = xb( 1 : nb );
    xb_y = xb( 1 + nb : 2*nb );
    
    %find x points within support
    xv_r = repmat( xv, [1,length(xb_x)]);
    xb_rx = repelem( xb_x, length(xv) );
    ind_supp_x = ( abs( xv_r - xb_rx ) <= supp * del );
    
    %Get index of these points for fluid grid
    i_supp = round( (xv_r(ind_supp_x) + offx)/del + 1/2 );
    
    %Get index of these points for IB:
    i_supp_xbx = round( ind_supp_x .* repelem( 1:nb, length(xv) ) );
    i_supp_xbx = i_supp_xbx( i_supp_xbx ~= 0 );
        
    %find y points within support
    yv_r = repmat( yv, [1,length(xb_y)]);
    xb_ry = repelem( xb_y, length(yv) );
    ind_supp_y = ( abs( yv_r - xb_ry ) <= supp * del );
    
    %Get y-index of these points for fluid grid
    j_supp = round( (yv_r(ind_supp_y) + offy) /del );
    
    %Get index of these points for IB:
    j_supp_xby = round( ind_supp_y .* repelem( 1:nb, length(yv) ) );
    j_supp_xby = j_supp_xby( j_supp_xby ~= 0 );
    
    %For each IB point, add nonzero weights...
    for j = 1 : nb
        
        %x-indices on IB corresponding to current body point
        ind_xbx = (i_supp_xbx == j);
        
        %x-indices on flow grid that are within support of IB point
        ind_x = i_supp( ind_xbx );
        
        %y-indices on IB corresponding to current body point
        ind_xby = (j_supp_xby == j);
        
        %y-indices on flow grid that are within support of IB point
        ind_y = j_supp( ind_xby );
        
        %Combine flow indices
        indvelx = repmat( ind_x, [1,length(ind_y)]);
        indvely = repelem( ind_y, length(ind_x) );
        
        %rows of ET to add to
        vely_ind = get_vely_ind( indvelx, indvely, 1, grid_parms );
        
        %columns to add to
        xb_ind = j * ones( size(vely_ind ) );
        
        %entries to put into columns
        del_h = delta_h( xv(indvelx), xb_x(xb_ind), del) .* ...
            delta_h( yv(indvely), xb_y(xb_ind), del) ;
        
        
        %Add to ET:
        %y index starts after all the x-vels
        n_add = get_velx_ind( m-1, n, 1, grid_parms );
        ET = ET + sparse( n_add + vely_ind, xb_ind + nb, ...
            del_h, nrows, ncols );

        
    end
    
%--




