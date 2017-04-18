function rhs = get_Lap_BCs( rhs, parms, mats )

%Add contributions to Laplacian from outer edge of domain (provided by
%coarser grid)

    if lev < mg
        %**Now for BCs

            %scaling factor to multiply gamma by:
            %   The 1/4 turns circ on coarser grid into circ on finer grid
            scl = 1/4;

            %Bottom part

                %indices for coarse grid corresponding to bottom part of fine
                %domain
                bottom = (n/4 - 1)*(m-1) + ( m/4 : 3*m/4 );

                %points that need to average coarser domain:
                rhs( 1:2:m-1 ) = rhs( 1:2:m-1 ) + ...
                    1/2 * scl * rhs( bottom(1:end-1), lev + 1) + ...
                    1/2 * scl * rhs( bottom(2:end), lev + 1);

                %points that don't need to average coarser domain:
                rhs( 2:2:m-2 ) = rhs( 2:2:m-2 ) + ...
                    scl * rhs( bottom(2:end-1), lev + 1);

            %Top part

                %indices on coarse grid corresponding to top part of fine grid
                top = (3*n/4 - 1)*(m-1) + ( m/4 : 3*m/4 );

                %indices for fine grid
                topf = (n-2)*(m-1) + (1 : m-1);

                %points that need to average coarser domain:
                Wgam( topf(1:2:end) ) = Wgam( topf(1:2:end) ) + ...
                    1/2 * scl * gamma( top(1:end-1), lev + 1) + ...
                    1/2 * scl * gamma( top(2:end), lev + 1);

                %points that don't need to average coarser domain:
                Wgam( topf(2:2:end-1) ) = Wgam( topf(2:2:end-1) ) + ...
                    scl * gamma( top(2:end-1), lev + 1);


            %left part

                %indices on coarse grid corresponding to left edge of fine grid
                left = (n/4 - 1)*(m-1) + ( m/4 : m-1 : (m-1)*n/2 + m/4 );

                %indices for fine grid
                leftf = 1 : m-1 : (m-1)*(n-2) + 1 ;

                %points that need to average coarser domain:
                Wgam( leftf(1:2:end) ) = Wgam( leftf(1:2:end) ) + ...
                    1/2 * scl * gamma( left(1:end-1), lev + 1) + ...
                    1/2 * scl * gamma( left(2:end), lev + 1);

                %points that don't need to average coarser domain:
                Wgam( leftf(2:2:end-1) ) = Wgam( leftf(2:2:end-1) ) + ...
                    scl * gamma( left(2:end-1), lev + 1);

            %right part

                %indices on coarse grid corresponding to right edge of fine grid
                right = (n/4 - 1)*(m-1) + ( 3*m/4 : m-1 : (m-1)*n/2 + 3*m/4 );

                %indices for fine grid
                rightf = m-1 : m-1 : (m-1)*(n-2) + m-1 ;

                %points that need to average coarser domain:
                Wgam( rightf(1:2:end) ) = Wgam( rightf(1:2:end) ) + ...
                    1/2 * scl * gamma( right(1:end-1), lev + 1) + ...
                    1/2 * scl * gamma( right(2:end), lev + 1);

                %points that don't need to average coarser domain:
                Wgam( rightf(2:2:end-1) ) = Wgam( rightf(2:2:end-1) ) + ...
                    scl * gamma( right(2:end-1), lev + 1);

            %!!
        %**
    end


