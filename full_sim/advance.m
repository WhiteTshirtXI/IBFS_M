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

        q = zeros( nq, parms.mg );
        gamma = zeros( ngam, parms.mg );
        
        %freestream flux that doesn't contribute to vorticity
        q0 = zeros( nq, parms.mg );
        for j = 1 : parms.mg

            %grid spacing on current grid
            hc = h * 2^(j-1);

            %write fluid velocity flux in body-fixed frame
            q0( :, j ) = -parms.U_body * hc;
        end
        soln.q0 = q0;
        
    else
        q = soln.q;
        gamma = soln.gamma;
        q0 = soln.q0;
        nonlin_prev = soln.nonlin_prev;
    end
%--
    

%--Loop through different gridlevels and get trial circ

    nonlin = zeros( ngam, parms.mg );
    for j = parms.mg : -1 : 1

        %**Build rhs (explicit part of time-stepping)

            %compute the nonlinear term 
            nonlin(:,j) = get_nonlin( gamma, q, q0, j, parms, mats );

            if it == 0
                nonlin_prev = nonlin;
            end

            %combine explicit Laplacian and nonlinear terms into a rhs
            rhs = (mats.I - dt/2 * mats.Lap) * gamma - ...
                3*dt/2 * nonlin + dt/2 * nonlin_prev;
            
            %add boundary conditions from Laplacian term
            rhs = get_Lap_BCs( rhs, parms, mats );

            %Store nonlinear solution for use in next time step
            soln.nonlin_prev = nonlin;
        %**

        %**trial circulation

            gamm_star = mats.invIdtLap( rhs ); 

        %**

    end
%--
    
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
        