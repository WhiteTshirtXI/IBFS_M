function [soln,parms,mats] = advance( it, parms, mats, soln )


%--define variables used frequently
    
    %grid spacing on finest grid
    h = parms.len / parms.m;
    %# of x-vel (flux) points
    nu = parms.n * (parms.m-1); 
    %# of y-vel (flux) points
    nv = parms.m * (parms.n-1);
    %Total # of vel (flux) points
    nq = nu + nv;
    %# of vort (circ) points
    ngam = (parms.m-1) * (parms.n-1);
    %# of streamfunction points (same as # of vort points)
    nstrm = ngam;
    %# of body points
    nb = parms.nb;
    %# of surface stress points
    nf = nb * 2;
    %time step size in physical time
    dt = parms.dt;
    %Grid spacing on IB
    ds = parms.ds;
    
%--

%--
%need to build and store matrix and its inverse if at first time step
    if ( it == parms.it_start & parms.stationary == 'T' )
        
        display('building and storing matrix for computing surface stresses.')

        mats.B = zeros( nf, nf );
        for j = 1 : nf

            z = zeros( nf, 1 );
            z(j) = 1;
            
            mats.B(:,j) = b_times( z, parms, mats );
            
        end

        display('Pre-computing inverse for surface stresses.')
        mats.Binv = mats.B \ eye( nf, nf ) ;
        display('done!!')

    end
    
%--

%--variables required to advance

    if (it == 0)

        q = zeros( nq, parms.mg );
        gamma = zeros( ngam, parms.mg );
        
        chi = zeros( 3 * nb, 1); soln.chi = chi;
        zeta = chi;
        zetad = zeta;
        fb_til_dt = zeros( 2*nb, 1);
        
        %freestream flux that doesn't contribute to vorticity
        q0 = zeros( nq, parms.mg );
        for j = 1 : parms.mg

            %grid spacing on current grid
            hc = h * 2^(j-1);

            %write fluid velocity flux in body-fixed frame
            q0( 1 : (parms.m-1)*parms.n, j ) = -parms.U_body * hc;
            
        end
        
        soln.q0 = q0;
        
    else
        q = soln.q;
        gamma = soln.gamma;
        q0 = soln.q0;
        nonlin_prev = soln.nonlin_prev;
        gamm_star = zeros( ngam, parms.mg );
        nonlin = zeros( ngam, parms.mg );
        rhs = nonlin;
        
        if parms.deform == 'T'
            chi = soln.chi;
            zeta = soln.zeta;
            zetad = soln.zetad;
            fb_til_dt =  soln.fb *  ds / h * dt;

        end
    end
%--
    

%--Loop through different gridlevels and get trial circ that doesn't
%  satisfy no slip BC on IB

    for j = parms.mg : -1 : 1
        
        %grid spacing on current grid
            hc = h * 2^(j-1);

        %**Build rhs (explicit part of time-stepping)

            %compute the nonlinear term 
            nonlin(:,j) = get_nonlin( gamma, q, q0, j, parms, mats );

            if it == 0
                nonlin_prev = nonlin;
            end
            
            %contribution of Laplacian term...
                rhsbc = zeros( ngam, 1 );
                
                if j < parms.mg
                    %from explicit treatment of circulation
                    rhsbc = get_Lap_BCs( rhsbc, gamma(:,j+1), j, parms );

                    %from current (trial) vorticity
                    rhsbc = get_Lap_BCs( rhsbc, gamm_star(:,j+1), j, parms );
                end
                
            %combine explicit Laplacian and nonlinear terms into a rhs
            rhs(:,j) = (mats.I - dt/2 * mats.Lap / (hc^2) ) * ...
                gamma(:,j) - 3*dt/2 * nonlin(:,j) + ...
                dt/2 * nonlin_prev(:,j) + dt/2 * rhsbc;
            
        %**
        
        
        
        %kick
        if it == 10
            rhs(24, 1) = rhs(24,1) + 1e-2;
        end
        
        
        
        
        
        
        
        

        %**trial circulation

            gamm_star(:,j) = Ainv( rhs(:,j), j, parms, mats ); 
            %(The j means we are evaluating Ainv at the first grid level)
            
        %**

    end
    
    %Store nonlinear solution for use in next time step
    soln.nonlin_prev = nonlin;
%--
    
%--Update body position and surface stresses 
   
    if parms.stationary == 'T'
        %Get surface stresses from fluid onto body 
        %NOTE: stresses are multiplied by dt and are off from the physical
        %      stress by a scaling factor of ds/h.

        %don't need all grid levels to get trial velocity on 1st grid level
        if parms.mg > 1
            q_star = circ2_st_vflx( gamm_star, 2, parms, mats );
        else
            q_star = circ2_st_vflx( gamm_star, 1, parms, mats );
        end

        fb_til_dt = mats.Binv * ( 1/h * mats.E * q_star(:,1) + ...
            1/h * mats.E * q0(:,1) ); 

        soln.fb = fb_til_dt * h / ds / dt ; %surface stress in physical units

        soln.CD( it + 1 ) = 2 * sum( soln.fb(1 : nb) ) * ds;
        soln.CL( it + 1 ) = 2 * sum( soln.fb(1 + nb : nf ) ) * ds;

    else %deforming body
        
        %--step 1: initialize
            err_fsi = 10.d0;
            tol_fsi = 1.d-5;

            %Update matrices that depend on flag position and velocity
            mats = get_mats( parms, mats, soln );

            %store old body position and velocity
            chiold = chi;
            zetaold = zeta;
            zetadold = zetad;

        %--
    
        
        while err_fsi > tol_fsi

            %--step 2: compute guess for surface stresses
                
                %Flag matrix that comes up a bunch
                mats.sol_mat = mats.K_flag + 4/(dt^2) * mats.M_flag;
                mats.sol_mat = mats.sol_mat \ eye( size( mats.sol_mat) );

                %--rhs contribution from constraint

                    r_cb = 2/dt * (chi - chiold) - zetaold;

                %--

                %--rhs contribution from structure equations

                    r_ub = mats.M_flag * (zetadold + 4/dt * zetaold + ...
                        4/(dt^2)* (chiold - chi) ) - mats.Fint;

                    r_ub = mats.sol_mat * r_ub;

                %--

                %--put rhs together and truncate to be size of surface stress
                %  vector

                    F_sm = -mats.Itilde_flag * ( 2/dt * r_ub + r_cb );
                    
                    
                    %Need to convert trial circ to vel flux
                    %(don't need all grid levs for trial vel on 1st grid
                    %lev)
                    if parms.mg > 1
                        q_star = circ2_st_vflx( gamm_star, 2, parms, mats );
                    else
                        q_star = circ2_st_vflx( gamm_star, 1, parms, mats );
                    end

                    rhsf = 1/h * mats.E * q_star(:,1) + ...
                        1/h * mats.E * q0(:,1) + F_sm;
                %--

                %--get surface stress
                
                    fb_til_dt = bicgstab( fb_til_dt, rhsf, parms, mats ) ;
                    %NOTE: stresses are multiplied by dt different from 
                    %       physical stress by a scaling factor of ds/h.
                    
                    soln.fb = fb_til_dt * h / ds / dt ; %surface stress in physical units

                %--

            %--

            %--step 3: compute increment in body position

                d_chi = r_ub + mats.sol_mat * mats.Q_flag * mats.W_flag ...
                    * fb_til_dt * ( h/ds ) / dt; 

                err_fsi = max(abs( d_chi ) ) / max(abs( chi ) );

            %--

            %--step 4: update body positions and velocities

                chi = chi + d_chi;
                zeta = -zetaold + 2/dt * (chi - chiold);
                zetad = 4/(dt^2) * ( chi - chiold) - 4/dt*zetaold - zetadold;

                soln.xb = parms.xb0 + (mats.Itilde_flag * chi)';
                
                %store flag quantities for next time step
                soln.chi = chi;
                soln.zeta = zeta;
                soln.zetad = zetad;

            %--

            %--step 5: update for next iteration if necessary

                if ( err_fsi >= tol_fsi )

                    %Update matrices that depend on flag position and vel
                    mats = get_mats( parms, mats, soln );

                end

            %--
       

        end
            
            soln.CD( it + 1 ) = 2 * sum( soln.fb(1 : nb) ) * ds;
            soln.CL( it + 1 ) = 2 * sum( soln.fb(1 + nb : nf ) ) * ds;
            
            if parms.inverted == 'T'
                soln.tip_disp( it+1 ) = soln.chi( 2 );
            else
                soln.tip_disp( it+1 ) = soln.chi( end - 1 );
            end
            
        %--
        
        
    end
    
   
%--

%--update circulation on fine grid to satisfy no-slip

    gamma = gamm_star;
    
    gamma(:,1) = gamm_star(:,1) - Ainv( mats.R * mats.ET ...
        * fb_til_dt, 1, parms, mats ) ;
    %Note we don't include BCs from coarse grid for Ainv because surface
    %stresses are compact
    
%--

%--Update circ on all grids

    for j = 2 : parms.mg

        gamma(:, j) = coarsify( gamma(:,j-1), gamma(:,j), parms );
        
    end
    soln.gamma = gamma;
    
%--

%--Get vel flux and streamfcn from circulation

    [soln.q, soln.s] = circ2_st_vflx( gamma, parms.mg, parms, mats );

%--

%--A few simulation quantities of interest
    
    %slip on IB
    soln.slip( it + 1 ) = 1/h * max( abs( mats.E * soln.q(:,1) + ...
        mats.E*q0(:,1) ) );
    
    %get cfl (u * dt / dx) :
    soln.cfl( it + 1 ) = max( abs( 1/(h^2) * soln.q(:,1) * dt ) ) ; 
    
%--
        
