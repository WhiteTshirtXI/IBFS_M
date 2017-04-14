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


% % %--
% % %need to build and store matrix and its inverse if at first time step
% %     if it == 0
% %         
% %         display('building and storing matrix inverse for computing surface stresses.')
% %         mats.B = mats.E * mats.M_vel * mats.C * mats.invRC( ...
% %             mats.invIdtLap( mats.R * mats.ET ) ) ;
% %         mats.Binv = mats.B \ eye( nf, nf ) ;
% %         display('done!!')
% % 
% %     end
% %     
% % %--

%--variables required to advance

    if (it == 0)

        q = zeros( nq, 1 );
        gamma = zeros( ngam, 1 );
        chi = zeros( 3 * nb, 1);
        zeta = chi;
        zetad = zeta;
        fb_til_dt = zeros( 2*nb, 1);
        
        soln.chi = chi;
        soln.fb = fb_til_dt;
        
        
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
        chi = soln.chi;
        zeta = soln.zeta;
        zetad = soln.zetad;
        fb_til_dt =  soln.fb *  ds / h * dt;
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




%add perturbation:
if it == 15
    rhs( 200 ) = rhs(200) + 1e-1;
end




%--trial circulation

    gamm_star = mats.invIdtLap( rhs ); 
    
%--

%--FSI loop

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
        
        
            %************NEED TO LOOK INTO WHETHER IT'S BETTER TO
            %PRE-COMPUTE AND STORE THE INVERSE TO THIS BAD BOY
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
               
                rhsf = mats.E * mats.M_vel * mats.C * ...
                ( mats.invRC( gamm_star ) ) + mats.E* mats.M_vel * q0 ...
                + F_sm;
            %--
            
            %--get surface stress
                
                fb_til_dt = bicgstab( fb_til_dt, rhsf, parms, mats ) ;
                %NOTE: stresses are multiplied by dt different from 
                %       physical stress by a scaling factor of ds/h.
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
            
        %--
            
        %--step 5: update for next iteration if necessary
        
            if (err_fsi >= tol_fsi )
               
                %Update matrices that depend on flag position and velocity
                mats = get_mats( parms, mats, soln );
                
            end
        
        %--
           
            
    end
    
    CD = 2 * sum( fb_til_dt(1 : nb) * h / ds /dt ) * ds
    CL = 2 * sum( fb_til_dt(nb + (1 : nb) ) * h / ds /dt ) * ds
    tip_disp = chi( 2 )
%--

%--update circulation to satisfy no-slip


    gamma = gamm_star - mats.invIdtLap( mats.R * mats.ET * fb_til_dt ) ;

%--
    
    soln.chi = chi;
    soln.zeta = zeta;
    soln.zetad = zetad;
    soln.q = mats.C * ( mats.invRC( gamma ) );
    soln.gamma = gamma;
    soln.fb = fb_til_dt * h / ds /dt ; %get surface stress in physical units
    soln.fb_rdst = mats.W_flag * soln.fb;
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
        