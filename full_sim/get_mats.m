function mats = get_mats( parms )

%Build and store matrices required for simulation

%---Curl matrices

    display('getting curl and rot matrices')    
    [mats.C, mats.R] = get_curls( parms ); 
    
    mats.RC = mats.R * mats.C;
    
% %     display('getting LU factorization of RC. Patience, grasshoppa...')
% %     [mats.LRC,mats.URC,mats.pRC,mats.qRC,mats.rRC] = lu(mats.RC);
% % 
% %     %get matrix function handle for solving linear system
% %     mats.invRC = @(x) mats.qRC*(mats.URC\(mats.LRC\(mats.pRC*(mats.rRC \ x ) ) ) );
   
%---

%---Scaling matrix that converts velocities to fluxes

    display('getting vel to flux scaling matrix')
    mats.M_vel = get_M_vel( parms );
    
    %also need its inverse:
    mats.M_vel_inv = mats.M_vel;
    mats.M_vel_inv( mats.M_vel_inv ~= 0 ) = ...
        1./mats.M_vel_inv( mats.M_vel_inv ~= 0 );
    
%---

%---Scaling matrix that converts vorticity to circulation

    display('getting vort to circ scaling matrix')
    mats.M_vort = get_M_vort( parms );
    
%---


%---Laplacian

    mats.Lap = 1/parms.Re * mats.RC ;
    
    %-Identity matrix used in conjunction with Laplacian

    mats.I = speye( size( mats.Lap ) );

    %-
    
    %eigenvalues of Laplacian
    ii = 1 : (m-1); jj = 1 : (n-1);
    [ii, jj] = meshgrid( ii, jj );
    mats.lam = 2 .* ( cos( pi .* ii ./ m ) + cos( pi .* jj ./ n ) - 2);
    
% %     display('getting LU factorization of (I+dt/2 * Lap). Patience, grasshoppa...')
% %     [mats.LLap,mats.ULap,mats.pLap,mats.qLap,mats.rLap] = lu( ...
% %         mats.I + parms.dt/2 * mats.Lap );
% % 
% %     %get matrix function handle for solving linear system
% %     mats.invIdtLap = @(x) mats.qLap*( mats.ULap\( mats.LLap\( ...
% %         mats.pLap*( mats.rLap \ x ) ) ) );
   

    
%---



%---matrices for nonlinear term

    display('getting Q (for nonlinear term)')
    mats.Q = get_Q( parms );
    
    mats.Q = mats.Q * mats.M_vel;
    
    display('getting W (for nonlinear term)')
    mats.W = get_W( parms );
   
    mats.W = mats.M_vel_inv * mats.W;
    
%---


%---Smearing and interpolation matrices

    display('getting smearing and interpolation matrices')
    parms.supp = 6; %currently using a delta function with 
                    %a support of 6 grid cells...
                    
    
    mats.ET = get_ET( parms );

    mats.E = mats.ET' ;
    
%---
    




