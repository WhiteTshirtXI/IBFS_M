function g = g_eval( parms, mats, soln )


%--Various variables 
    %# of x-vel (flux) points
    nu = get_velx_ind( parms.m-1, parms.n, parms.mg, parms );
    %# of y-vel (flux) points
    nv = get_vely_ind( parms.m, parms.n-1, parms.mg, parms );
    %Total # of vel (flux) points
    nq = nu + nv;
    %# of vort (circ) points
    ngam = get_vort_ind( parms.m, parms.n-1, parms.mg, parms );
    %# of surface stress points
    nf = 2*parms.nb;
    %grid spacing on finest grid (used a lot)
    h = parms.len / parms.m;
%--

%--specify matrices for ease in ensuing code

    C = mats.C; R = mats.R; M_vel = mats.M_vel; Lap = mats.Lap;
    Q = mats.Q; W = mats.W; ET = mats.ET; E = mats.E;
    I = mats.I; RC = mats.RC; I_til = mats.I_til;
    I_til_s2vort = mats.I_til_s2vort;

%--

%--soln vars

    s = soln.s; fb = soln.fb; q0 = parms.q0;

%--

%--build g
    g = zeros( ngam + nf, 1);
    
    
    %rows corresponding to momentum equations
    diagWgam = sparse( 1:nq, 1:nq, W * R * C * s, nq, nq ); 
    
    g( 1 : ngam ) = -Lap * RC * s - R * ( diagWgam * ( Q* (C*s + q0) ) ) ...
             - R * ET * fb;
        
    %modify to account for Neumann BCs:
    g( 1 : ngam ) = I_til * g(1 : ngam) + I_til_s2vort*(RC*s);
    
    g( ngam + 1 : ngam + nf ) = E * M_vel * (C*s + q0 ) ;
%--






