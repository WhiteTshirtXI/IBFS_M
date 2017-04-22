function nonlin_v = get_nonlin( gamma, q, q0, lev, parms, mats )

%Build the nonlinear term. Without accounting for BCs, the answer is 
%nonlin = mats.R * ( ( mats.W * gamma ) .* ( mats.Q * (q + q0) ) );

%However, we need to modify the W*gamma term and the Q*(q + q0) term to
%account for circulation and velocity flux terms from the coarse grid on
%the fine grid. (Note that R does not need to be modified).


%--parameters used in this function

    m = parms.m; n = parms.n; mg = parms.mg;
    hf = parms.len / m; %grid spacing on finest level
    
    hc = hf * 2^( lev - 1); %grid spacing on current level
 
%--

%--Build Wgamma and Q(q + q0) without accounting for BCs

    %the 1/hc^2 term is to convert circ to vort
    Wgam = 1/(hc^2) * mats.W * gamma(:, lev );

    qq0 = q + q0;
    Qqq0 = mats.Q*qq0(:, lev);
    
%--



%--Now take BCs into account:
    
    if lev < mg
        
        %** Wgamma term
        
            %scaling factor to multiply gamma by:
            %   The 1/4 turns circ on coarser grid into circ on finer grid
            %   The 1/2 is because this term is part of an average
            %   The hc^2 term converts circ to vort
            scl = 1/2 * 1/4 / hc^2;

            %!!x-velocity block (bottom and top part of domain):

            %Bottom part

                %indices for coarse grid corresponding to bottom part of fine
                %domain
                bottom = (n/4 - 1)*(m-1) + ( m/4 : 3*m/4 );

                %points that need to average coarser domain:
                Wgam( 1:2:m-1 ) = Wgam( 1:2:m-1 ) + ...
                    1/2 * scl * gamma( bottom(1:end-1), lev + 1) + ...
                    1/2 * scl * gamma( bottom(2:end), lev + 1);

                %points that don't need to average coarser domain:
                Wgam( 2:2:m-2 ) = Wgam( 2:2:m-2 ) + ...
                    scl * gamma( bottom(2:end-1), lev + 1);

            %Top part

                %indices on coarse grid corresponding to top part of fine grid
                top = (3*n/4 - 1)*(m-1) + ( m/4 : 3*m/4 );

                %indices for fine grid
                topf = (n-1)*(m-1) + (1 : m-1);

                %points that need to average coarser domain:
                Wgam( topf(1:2:end) ) = Wgam( topf(1:2:end) ) + ...
                    1/2 * scl * gamma( top(1:end-1), lev + 1) + ...
                    1/2 * scl * gamma( top(2:end), lev + 1);

                %points that don't need to average coarser domain:
                Wgam( topf(2:2:end-1) ) = Wgam( topf(2:2:end-1) ) + ...
                    scl * gamma( top(2:end-1), lev + 1);


            %!!    

            %!!y-velocity block (contributions from left and right edges)
            
            %indices for y-vel start after x velocities
            nadd = (m-1)*n;
            
            %left part

                %indices on coarse grid corresponding to left edge of fine grid
                left = (n/4 - 1)*(m-1) + ( m/4 : m-1 : (m-1)*n/2 + m/4 );

                %indices for fine grid
                leftf = nadd + ( 1 : m : m*(n-2) + 1 );

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
                rightf = nadd + ( m : m : m*(n-2) + m);

                %points that need to average coarser domain:
                Wgam( rightf(1:2:end) ) = Wgam( rightf(1:2:end) ) + ...
                    1/2 * scl * gamma( right(1:end-1), lev + 1) + ...
                    1/2 * scl * gamma( right(2:end), lev + 1);

                %points that don't need to average coarser domain:
                Wgam( rightf(2:2:end-1) ) = Wgam( rightf(2:2:end-1) ) + ...
                    scl * gamma( right(2:end-1), lev + 1);

            %!!
        %**
        
        %** Qqq0 term
        
            %scaling factor to multiply qq0 by:
            %   The 1/2 converts the velocity flux on the coarser grid to a
            %   flux on the finer grid
            scl = 1/2 ;

            %!!x-velocity block (contributions to top and bottom edges)
            %   averages y-velocities
            
            %Bottom part

                %indices for coarse grid corresponding to bottom part of fine
                %domain
                bottom = nadd + (n/4 - 1) * m + ( m/4 + 1 : 3*m/4 );

                %points that need to average coarser domain:
                Qqq0( 2:2:m-2 ) = Qqq0( 2:2:m-2 ) - ...
                    1/4 * scl * qq0( bottom(1:end-1), lev + 1) - ...
                    1/4 * scl * qq0( bottom(2:end), lev + 1);

                %points that don't need to average coarser domain:
                Qqq0( 1:2:m-1 ) = Qqq0( 1:2:m-1 ) - ...
                    1/2 * scl * qq0( bottom, lev + 1);

            %Top part

                %indices on coarse grid corresponding to top part of fine grid
                top = nadd + (3*n/4 - 1)*m + ( m/4 + 1 : 3*m/4 );

                %indices for fine grid
                topf = (n-1)*(m-1) + (1 : m-1);

                %points that need to average coarser domain:
                Qqq0( topf(2:2:end-1) ) = Qqq0( topf(2:2:end-1) ) - ...
                    1/4 * scl * qq0( top(1:end-1), lev + 1) - ...
                    1/4 * scl * qq0( top(2:end), lev + 1);

                %points that don't need to average coarser domain:
                Qqq0( topf(1:2:end) ) = Qqq0( topf(1:2:end) ) - ...
                    1/2 * scl * qq0( top, lev + 1);


            %!!    

            %!!y-velocity block (averages x-velocities)
            
            %left part

                %indices on coarse grid corresponding to left edge of fine grid
                left = ( n/4 )*(m-1) + ( m/4 : m-1 : (m-1)*(n/2-1) + m/4 );

                %indices for fine grid
                leftf = nadd + (1 : m : m*(n-2) + 1 );

                %points that need to average coarser domain:
                Qqq0( leftf(2:2:end-1) ) = Qqq0( leftf(2:2:end-1) ) + ...
                    1/4 * scl * qq0( left(1:end-1), lev + 1) + ...
                    1/4 * scl * qq0( left(2:end), lev + 1);

                %points that don't need to average coarser domain:
                Qqq0( leftf(1:2:end) ) = Qqq0( leftf(1:2:end) ) + ...
                    1/2 * scl * qq0( left, lev + 1);

            %right part

                %indices on coarse grid corresponding to right edge of fine grid
                right = ( n/4 )*(m-1) + ( 3*m/4 : m-1 : (m-1)*(n/2-1) + 3*m/4 );

                %indices for fine grid
                rightf = nadd + ( m : m : m*(n-2) + m ) ;

                %points that need to average coarser domain:
                Qqq0( rightf(2:2:end-1) ) = Qqq0( rightf(2:2:end-1) ) + ...
                    1/4 * scl * qq0( right(1:end-1), lev + 1) + ...
                    1/4 * scl * qq0( right(2:end), lev + 1);

                %points that don't need to average coarser domain:
                Qqq0( rightf(1:2:end) ) = Qqq0( rightf(1:2:end) ) + ...
                    1/2 * scl * qq0( right, lev + 1);

            %!!
        %**
    end
%--




 nonlin_v = mats.R * ( Wgam .*  Qqq0  );



