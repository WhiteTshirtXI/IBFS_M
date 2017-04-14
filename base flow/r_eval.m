function r = r_eval( parms, mats, soln )


%--Various variables 
    %# of x-vel (flux) points
    nu = get_velx_ind( parms.m-1, parms.n, parms.mg, parms );
    %# of y-vel (flux) points
    nv = get_vely_ind( parms.m, parms.n-1, parms.mg, parms );
    %Total # of vel (flux) points
    nq = nu + nv;
    %# of vort (circ) points
    ngam = get_vort_ind( parms.m-1, parms.n-1, parms.mg, parms );
    %# of surface stress points
    nb = parms.nb;
    nf = 2*parms.nb;
    %grid spacing on finest grid (used a lot)
    h = parms.len / parms.m;
%--

%--specify matrices for ease in ensuing code

    C = mats.C; R = mats.R; M_vel = mats.M_vel; Lap = mats.Lap;
    Q = mats.Q; W = mats.W; ET = mats.ET; E = mats.E;
    I = mats.I; RC = mats.RC; M_fl = mats.M_flag; Q_fl = mats.Q_flag;
    W_fl = mats.W_flag; Itil_fl = mats.Itilde_flag; K_fl = mats.K_flag;
    J_Ef = mats.J_Ef; J_Wf = mats.J_Wf; J_Es = mats.J_Es;
    Fint = mats.Fint;

%--

%--soln vars

    s = soln.s; fb = soln.fb; q0 = parms.q0;
    chi = soln.chi; zeta = soln.zeta;
    
%--

%--build r
    r = zeros( ngam + 6*nb + nf, 1);
    
    
    %rows corresponding to fluid momentum equations
    diagWgam = sparse( 1:nq, 1:nq, W * R * C * s, nq, nq ); 
    
    r( 1 : ngam ) = -Lap * RC * s - R * ...
        ( diagWgam * ( Q* (C*s + q0) ) ) - R * ET * fb;
    
    %rows for structural momentum
%     r( ngam + (1 : 3*nb ) ) = -K_fl*chi + Q_fl * W_fl * fb;
    r( ngam + (1 : 3*nb ) ) = -Fint + Q_fl * W_fl * fb;
    
%     show = [-K_fl*chi -Fint Q_fl * W_fl * fb]
%     
%     pause
    
    %rows for body position
    r( ngam + 3*nb + (1 : 3*nb ) ) = zeta;
    
    %rows for no-slip BC
    r( ngam + 6*nb + (1 : nf) ) = E * M_vel * (C*s + q0 ) - ...
        Itil_fl* zeta ;
%--






