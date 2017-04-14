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
    ngam = get_vort_ind( parms.m, parms.n-1, parms.mg, parms );
    %# of body points
    nb = parms.nb;
    %# of surface stress points
    nf = nb * 2;
%--

%--specify matrices for ease in ensuing code

    C = mats.C; R = mats.R; M_vel = mats.M_vel; Lap = mats.Lap;
    Q = mats.Q; W = mats.W; ET = mats.ET; E = mats.E;
    I = mats.I; RC = mats.RC; I_til = mats.I_til; 
    I_til_s2vort = mats.I_til_s2vort;

%--

mats.Dg = sparse( ngam + nf, ngam + nf ); 

%--(1,1) block

    diagWgam = sparse( 1:nq, 1:nq, W * R * C * soln.s, nq, nq ); 
    diagQq = sparse( 1:nq, 1:nq, Q * ( C * soln.s + parms.q0), nq, nq );
    
    B11 = -Lap * RC - R * ( diagWgam * (Q * C) + diagQq * W * R * C );
    
    %modify block to account for Neumann BCs:
    B11 = I_til * B11 + I_til_s2vort * RC;
    
%     B11 = I_til * B11 + ( I - I_til ) * RC;

    
%     max(max(abs( I_til_s2vort - ( I - I_til ) ) ) )
    
    [i,j,s] = find( B11 );
    
    mats.Dg = mats.Dg + sparse( i,j,s, ngam+nf, ngam+nf );

%--

%--(1,2) block

    B12 = -R * ET;

    %modify block to account for Neumann BCs:
    B12 = I_til * B12;
    
    [i,j,s] = find( B12 );
    
    mats.Dg = mats.Dg + sparse( i, j + ngam , s , ngam+nf, ngam+nf );
%--

%--(2,1) block

    B21 = E * M_vel * C;

    [i,j,s] = find( B21 );
    
    mats.Dg = mats.Dg + sparse( i + ngam, j , s , ngam+nf, ngam+nf );
%--

JAC_time = toc;

% display('        done.')
% display(['        done. Spent ',num2str(JAC_time),' secs.'])
