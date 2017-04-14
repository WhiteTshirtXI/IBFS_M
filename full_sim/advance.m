function [soln,parms,mats] = advance( it, parms, mats, soln )


%--define variables used frequently
    
    %grid spacing on finest grid
    h = parms.len / parms.m;
    %# of x-vel (flux) points
    nu = get_velx_ind( parms.m-1, parms.n, parms.mg, parms );
    %# of y-vel (flux) points
    nv = get_vely_ind( parms.m, parms.n-1, parms.mg, parms );
    %Total # of vel (flux) points
    nq = nu + nv;
    %# of vort (circ) points
    ngam = get_vort_ind( parms.m-1, parms.n-1, parms.mg, parms );
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
    if it == 0
        
        display('building and storing matrix inverse for computing surface stresses.')
        mats.B = mats.E * mats.M_vel * mats.C * mats.invRC( ...
            mats.invIdtLap( mats.R * mats.ET ) ) ;
        mats.Binv = mats.B \ eye( nf, nf ) ;
        display('done!!')

    end
    
%--

%--variables required to advance

    if (it == 0)

        q = zeros( nq, 1 );
        gamma = zeros( ngam, 1 );
        
        %freestream flux that doesn't contribute to vorticity
        q0 = zeros( nq, 1 );
        for j = 1 : parms.mg

            %grid spacing on current grid
            hc = h * 2^(j-1);

            %index of x-vel flux points for current grid
            if j == 1
                ind_x = 1 : get_velx_ind( parms.m-1, parms.n, j, parms );
            else
                ind_x = 1 + get_velx_ind( parms.m-1, parms.n, j-1, parms ) : ...
                    get_velx_ind( parms.m-1, parms.n, j, parms );
            end

            %write fluid velocity flux in body-fixed frame
            q0( ind_x ) = -parms.U_body * hc;
        end
        soln.q0 = q0;
        
    else
        q = soln.q;
        gamma = soln.gamma;
        q0 = soln.q0;
        nonlin_prev = soln.nonlin_prev;
    end
%--
    
   
%--Build rhs (explicit part of time-stepping)

    %compute the nonlinear term 
    %           R * ( Q*q ).* (W * gamma) 
           
    nonlin = mats.R * ( ( mats.W * gamma ) .* ( mats.Q * (q + q0) ) );
    
    if it == 0
        nonlin_prev = nonlin;
    end

    %combine explicit Laplacian and nonlinear terms into a rhs
    rhs = (mats.I - dt/2 * mats.Lap) * gamma - ...
        3*dt/2 * nonlin + dt/2 * nonlin_prev;
        
    
    %Store nonlinear solution for use in next time step
    soln.nonlin_prev = nonlin;
%--

%--trial circulation

    gamm_star = mats.invIdtLap( rhs ); 
    
%--compute surface stress that corrects part of trial circulation not
%  satisfying the no-slip BCs
    
    %Get surface stresses from fluid onto body 
    %NOTE: stresses are multiplied by dt and are off from the physical
    %      stress by a scaling factor of ds/h.
    
    fb_til_dt = mats.Binv * ( mats.E * mats.M_vel * mats.C * ...
        ( mats.invRC( gamm_star ) ) + mats.E* mats.M_vel * q0 ); 
        %the q0 contribution is the part that gets 
        %removed by taking the curl.

%--

%--update circulation to satisfy no-slip


    gamma = gamm_star - mats.invIdtLap( mats.R * mats.ET * fb_til_dt ) ;

%--
    
    
    soln.q = mats.C * ( mats.invRC( gamma ) );
    soln.gamma = gamma;
    soln.fb = fb_til_dt * h / ds /dt ; %get surface stress in physical units
    soln.CD( it + 1 ) = 2 * sum( soln.fb(1 : nb) ) * h;
    soln.CL( it + 1 ) = 2 * sum( soln.fb(1 + nb : nf ) ) * h;
    soln.slip( it + 1 ) = max( abs( mats.E * mats.M_vel * soln.q - ...
        mats.E*mats.M_vel*q0 ) );
    
    %get cfl:
    soln.cfl( it + 1 ) = max( abs( mats.M_vel.^2 * soln.q * dt ) ) ; %udt / dx
    
    if mod( it, 100) == 0
        cfl = soln.cfl( it + 1 )
        CD = soln.CD( it + 1 )
    end
        