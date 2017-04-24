function mats = get_mats_preproc( parms, soln )

%Build and store matrices required for simulation

m = parms.m; n = parms.n;


%---Curl matrices

    display('getting curl and rot matrices')    
    [mats.C, mats.R] = get_curls( parms ); 
    
    mats.RC = mats.R * mats.C;
    
    %eigenvalues of RC (negative of the evals of the 5point stencil Lap)
    ii = 1 : (m-1); jj = 1 : (n-1);
    [ii, jj] = meshgrid( ii, jj );
    mats.lam = -2 .* ( cos( pi .* ii ./ m ) + cos( pi .* jj ./ n ) - 2);

%---


%---Laplacian

    mats.Lap = 1/parms.Re * mats.RC ;
    
    %-Identity matrix used in conjunction with Laplacian

    mats.I = speye( size( mats.Lap ) );

    %-
    
%---



%---matrices for nonlinear term

    display('getting Q (for nonlinear term)')
    mats.Q = get_Q( parms );
    
    display('getting W (for nonlinear term)')
    mats.W = get_W( parms );
   
    
%---

%---Smearing and interpolation matrices

    %only build and store for stationary bodies. Otherwise they will be
    %built in advance.m
    if parms.stationary == 'T'
        display('getting smearing and interpolation matrices')
        parms.supp = 6; %currently using a delta function with 
                        %a support of 6 grid cells...


        mats.ET = get_ET( soln.xb, parms );

        mats.E = mats.ET' ;
        
    end

%---
    



