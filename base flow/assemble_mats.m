function mats = assemble_mats( parms, mats, soln )

%Assemble Jacobian matrix for full system using natrices obtained from
%get_mats.


display('      Updating Jacobian...')

tic
%--Various variables 
    %# of x-vel (flux) points
    nu = get_velx_ind( parms.m-1, parms.n, parms.mg, parms );
    %# of y-vel (flux) points
    nv = get_vely_ind( parms.m, parms.n-1, parms.mg, parms );
    %Total # of vel (flux) points
    nq = nu + nv;
    %# of vort (circ) points
    ngam = get_vort_ind( parms.m-1, parms.n-1, parms.mg, parms );
    %# of body points
    nb = parms.nb;
    %# of surface stress points
    nf = nb * 2;
%--

%--specify matrices for ease in ensuing code

    C = mats.C; R = mats.R; M_vel = mats.M_vel; Lap = mats.Lap;
    Q = mats.Q; W = mats.W; ET = mats.ET; E = mats.E;
    I = mats.I; RC = mats.RC; M_fl = mats.M_flag; Q_fl = mats.Q_flag;
    W_fl = mats.W_flag; Itil_fl = mats.Itilde_flag; K_fl = mats.K_flag;
    J_Ef = mats.J_Ef; J_Wf = mats.J_Wf; J_Es = mats.J_Es;

%--

ndof = ngam + 6*nb + nf;
mats.Dr = sparse( ndof, ndof ); 

%--(1,1) block

    diagWgam = sparse( 1:nq, 1:nq, W * R * C * soln.s, nq, nq ); 
    diagQq = sparse( 1:nq, 1:nq, Q * ( C * soln.s + parms.q0), nq, nq );
    
    B11 = -Lap * RC - R * ( diagWgam * (Q * C) + diagQq * W * RC );
    
    [i,j,s] = find( B11 );
    
    mats.Dr = mats.Dr + sparse( i,j,s, ndof, ndof );

%--

%--(1,3) block

    B13 = -J_Ef;
    
    [i,j,s] = find( B13 );
    
    mats.Dr = mats.Dr + sparse( i,j + ngam + 3*nb ,s, ndof, ndof );

%--

%--(1,4) block

    B14 = -R * ET;

    [i,j,s] = find( B14 );
    
    mats.Dr = mats.Dr + sparse( i, j + ngam + 6*nb , s , ndof, ndof );
%--

%--(2,3) block

    B23 = -K_fl + J_Wf;
    
    [i,j,s] = find( B23 );
    
    mats.Dr = mats.Dr + sparse( i + ngam, j + ngam + 3*nb , s , ndof, ndof );
%--

%--(2,4) block

    B24 = Q_fl * W_fl;
    
    [i,j,s] = find( B24 );
    
    mats.Dr = mats.Dr + sparse( i + ngam, j + ngam + 6*nb , s , ndof, ndof );
%--

%--(3,2) block

    B32 = speye( 3*nb );
    
    [i,j,s] = find( B32 );
    
    mats.Dr = mats.Dr + sparse( i + ngam + 3*nb, j + ngam, s, ndof, ndof );
%--

%--(4,1) block

    B41 = E * M_vel * C;

    [i,j,s] = find( B41 );
    
    mats.Dr = mats.Dr + sparse( i + ngam + 6*nb, j , s , ndof, ndof );
%--

%--(4,2) block

    B42 = -Itil_fl;
    
    [i,j,s] = find( B42 );
    
    mats.Dr = mats.Dr + sparse( i + ngam + 6*nb, j + ngam , s , ndof, ndof );

%--

%--(4,3) block

    B43 = J_Es;

    [i,j,s] = find( B43 );
    
    mats.Dr = mats.Dr + sparse( i + ngam + 6*nb, j + ngam + 3*nb, s , ndof, ndof );

%--

JAC_time = toc;

% display('        done.')
% display(['        done. Spent ',num2str(JAC_time),' secs.'])
