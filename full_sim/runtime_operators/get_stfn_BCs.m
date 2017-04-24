function rhsbc = get_stfn_BCs( rhsbc, rhs, parms )

%BCs used for RCinv (taken from the streamfunction on the coarser grid)


m = parms.m; n = parms.n; mg = parms.mg;

%Add contributions to Laplacian from outer edge of domain (provided by
%coarser grid)

%**Now for BCs

    %Bottom part

        %indices for coarse grid corresponding to bottom part of fine
        %domain
        bottom = (n/4 - 1)*(m-1) + ( m/4 : 3*m/4 );

        %points that need to average coarser domain:
        rhsbc.bott( 1:2:m-1 ) = rhsbc.bott( 1:2:m-1 ) + ...
            1/2  * rhs( bottom(1:end-1) ) + ...
            1/2  * rhs( bottom(2:end) );

        %points that don't need to average coarser domain:
        rhsbc.bott( 2:2:m-2 ) = rhsbc.bott( 2:2:m-2 ) + ...
             rhs( bottom(2:end-1) );

    %Top part

        %indices on coarse grid corresponding to top part of fine grid
        top = (3*n/4 - 1)*(m-1) + ( m/4 : 3*m/4 );

        %indices for fine grid
        topf = (n-2)*(m-1) + (1 : m-1);

        %points that need to average coarser domain:
        rhsbc.top( topf(1:2:end) ) = rhsbc.top( topf(1:2:end) ) + ...
            1/2  * rhs( top(1:end-1) ) + ...
            1/2  * rhs( top(2:end) );

        %points that don't need to average coarser domain:
        rhsbc.top( topf(2:2:end-1) ) = rhsbc.top( topf(2:2:end-1) ) + ...
             rhs( top(2:end-1) );


    %left part

        %indices on coarse grid corresponding to left edge of fine grid
        left = (n/4 - 1)*(m-1) + ( m/4 : m-1 : (m-1)*n/2 + m/4 );

        %indices for fine grid
        leftf = 1 : m-1 : (m-1)*(n-2) + 1 ;

        %points that need to average coarser domain:
        rhsbc.left( leftf(1:2:end) ) = rhsbc.left( leftf(1:2:end) ) + ...
            1/2  * rhs( left(1:end-1) ) + ...
            1/2  * rhs( left(2:end) );

        %points that don't need to average coarser domain:
        rhsbc.left( leftf(2:2:end-1) ) = rhsbc.left( leftf(2:2:end-1) ) + ...
             rhs( left(2:end-1) );

    %right part

        %indices on coarse grid corresponding to right edge of fine grid
        right = (n/4 - 1)*(m-1) + ( 3*m/4 : m-1 : (m-1)*n/2 + 3*m/4 );

        %indices for fine grid
        rightf = m-1 : m-1 : (m-1)*(n-2) + m-1 ;

        %points that need to average coarser domain:
        rhsbc.right( rightf(1:2:end) ) = rhsbc.right( rightf(1:2:end) ) + ...
            1/2  * rhs( right(1:end-1) ) + ...
            1/2  * rhs( right(2:end) );

        %points that don't need to average coarser domain:
        rhsbc.right( rightf(2:2:end-1) ) = rhsbc.right( rightf(2:2:end-1) ) + ...
             rhs( right(2:end-1) );

    %!!
%**



