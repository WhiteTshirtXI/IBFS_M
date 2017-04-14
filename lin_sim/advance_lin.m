function [soln_lin, mats] = advance_lin( it, soln_lin, parms, mats, soln )


%grid spacing on finest grid (used a lot)
h = parms.len / parms.m;

%--define number of points for various quantities
    
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

%--specify matrices for ease in ensuing code

    C = mats.C; R = mats.R; M_vel = mats.M_vel; Lap = mats.Lap;
    Q = mats.Q; W = mats.W; ET = mats.ET; E = mats.E;
    I = mats.I; RC = mats.RC;

%--

%--
%need to build and store matrix and its inverse if at first time step
    if it == 0
        
        display('building and storing matrix inverse for computing surface stresses.')
        B1 = ( I + dt/2 * Lap ) \ ( R * ET );
        B2 = RC \ B1;
        mats.B = E * M_vel * C * B2;
        mats.Binv = mats.B \ eye( nf, nf ) ;
        display('done!!')

    end
    Binv = mats.Binv; 
    
%--

%--variables required to advance

    %base variables
    qbase = soln.q;
    gammabase = soln.gamma;
    q0base = parms.q0;

    if (it == 0)

        q = zeros( nq, 1 );
        gamma = zeros( ngam, 1 );
        soln_lin.q0 = zeros( nq, 1);
        q0 = soln_lin.q0;
        
    else
        q = soln_lin.q;
        gamma = soln_lin.gamma;
        q0 = soln_lin.q0;
        nonlin_prev = soln_lin.nonlin_prev;
    end
%--
    
   
%--Build rhs (explicit part of time-stepping)

    %compute the nonlinear term 
    %           R * (W * gamma) .* ( Q*q )
           
    nonlin = R * ( ( W * gammabase ) .* ( Q * (q + q0) ) ) + ...
        R * ( ( W * gamma ) .* ( Q * (qbase + q0base) ) );
    
    if it == 0
        nonlin_prev = nonlin;
    end

    %combine explicit Laplacian and nonlinear terms into a rhs
    rhs = (I - dt/2 * Lap) * gamma - ...
        3*dt/2 * nonlin + dt/2 * nonlin_prev;
    
    %kick the nonlinear solver
    if it == 3
        rhsbf = zeros( ngam, 1 );
        rhsbf( get_vort_ind( parms.m/2, parms.n/4+1 : parms.n/2 , 1, parms ) ) = ...
            h;
        rhs = rhs + rhsbf;
    end
        
    
    %Store nonlinear solution for use in next time step
    soln_lin.nonlin_prev = nonlin;
%--

%--trial circulation

    gamm_star = ( I + dt/2 * Lap ) \ rhs;
    
%--compute surface stress that corrects part of trial circulation not
%  satisfying the no-slip BCs
    
    %Get surface stresses from fluid onto body 
    %NOTE: stresses are multiplied by dt and are off from the physical
    %      stress by a scaling factor of ds/h.
    
    fb_til_dt = Binv * ( E * M_vel * C * (RC \ gamm_star) + ...
        E * M_vel * q0 ); %the q0 contribution is the part that gets 
                           %removed by taking the curl.

%--

%--update circulation to satisfy no-slip

    gamma = gamm_star - ( ( I + dt/2 *Lap ) \ ( R * ET * fb_til_dt ) );

%--
    
    
    
    soln_lin.q = C * ( RC \ gamma );
    soln_lin.gamma = gamma;
    soln_lin.fb = fb_til_dt * h / ds /dt ; %get surface stress in physical units
    soln_lin.CD( it + 1 ) = 2 * sum( soln_lin.fb(1 : nb) ) * h;
    soln_lin.CL( it + 1 ) = 2 * sum( soln_lin.fb(1 + nb : nf ) ) * h;
    soln_lin.slip( it + 1 ) = max( abs( E * M_vel * soln_lin.q - E*M_vel*q0 ) );
    
    %get cfl:
    soln_lin.cfl( it + 1 ) = max( abs( M_vel.^2 * soln_lin.q * dt ) ) ; %udt / dx
    
    if mod( it, 10) == 0
        cfl = soln_lin.cfl( it + 1 )
        CD = soln_lin.CD( it + 1 )
    end
        