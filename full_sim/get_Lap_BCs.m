function rhsbc = get_Lap_BCs( rhsbc, rhs, lev, parms )

%BCs used for Ainv


m = parms.m; n = parms.n; 
h = parms.len/ parms.m;

%grid spacing on current grid
hc = h * 2^(lev-1);


%Add contributions to Laplacian from outer edge of domain (provided by
%coarser grid)

%**Now for BCs

    %scaling factor to multiply gamma by:
    %   The 1/4 turns circ on coarser grid into circ on finer grid
    %   The 1/hc^2 is because the Laplacian is a second deriv operator
    %   The 1/Re number is the physical scaling for the viscous Lap 
    scl = 1/4 / (hc^2) / parms.Re;

    %Bottom part

        %indices for coarse grid corresponding to bottom part of fine
        %domain
        bottom = (n/4 - 1)*(m-1) + ( m/4 : 3*m/4 );

        %points that need to average coarser domain:
        rhsbc( 1:2:m-1 ) = rhsbc( 1:2:m-1 ) + ...
            1/2 * scl * rhs( bottom(1:end-1)) + ...
            1/2 * scl * rhs( bottom(2:end));

        %points that don't need to average coarser domain:
        rhsbc( 2:2:m-2 ) = rhsbc( 2:2:m-2 ) + ...
            scl * rhs( bottom(2:end-1));

    %Top part

        %indices on coarse grid corresponding to top part of fine grid
        top = (3*n/4 - 1)*(m-1) + ( m/4 : 3*m/4 );

        %indices for fine grid
        topf = (n-2)*(m-1) + (1 : m-1);

        %points that need to average coarser domain:
        rhsbc( topf(1:2:end) ) = rhsbc( topf(1:2:end) ) + ...
            1/2 * scl * rhs( top(1:end-1)) + ...
            1/2 * scl * rhs( top(2:end));

        %points that don't need to average coarser domain:
        rhsbc( topf(2:2:end-1) ) = rhsbc( topf(2:2:end-1) ) + ...
            scl * rhs( top(2:end-1));


    %left part

        %indices on coarse grid corresponding to left edge of fine grid
        left = (n/4 - 1)*(m-1) + ( m/4 : m-1 : (m-1)*n/2 + m/4 );

        %indices for fine grid
        leftf = 1 : m-1 : (m-1)*(n-2) + 1 ;

        %points that need to average coarser domain:
        rhsbc( leftf(1:2:end) ) = rhsbc( leftf(1:2:end) ) + ...
            1/2 * scl * rhs( left(1:end-1)) + ...
            1/2 * scl * rhs( left(2:end));

        %points that don't need to average coarser domain:
        rhsbc( leftf(2:2:end-1) ) = rhsbc( leftf(2:2:end-1) ) + ...
            scl * rhs( left(2:end-1));

    %right part

        %indices on coarse grid corresponding to right edge of fine grid
        right = (n/4 - 1)*(m-1) + ( 3*m/4 : m-1 : (m-1)*n/2 + 3*m/4 );

        %indices for fine grid
        rightf = m-1 : m-1 : (m-1)*(n-2) + m-1 ;

        %points that need to average coarser domain:
        rhsbc( rightf(1:2:end) ) = rhsbc( rightf(1:2:end) ) + ...
            1/2 * scl * rhs( right(1:end-1)) + ...
            1/2 * scl * rhs( right(2:end));

        %points that don't need to average coarser domain:
        rhsbc( rightf(2:2:end-1) ) = rhsbc( rightf(2:2:end-1) ) + ...
            scl * rhs( right(2:end-1));

    
%**



