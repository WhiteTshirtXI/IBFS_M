function mats = get_mats( parms )

%Build and store matrices required for simulation

%---Curl matrices

    display('getting curl and rot matrices')    
    [mats.C, mats.R] = get_curls( parms ); 
    
    mats.RC = mats.R * mats.C;
   
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
    mats.Lap = 1/parms.Re * mats.R * mats.C * mats.M_vort ;
    
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
    

%---Identity matrix used in conjunction with Laplacian

    mats.I = speye( size( mats.Lap ) );

%---

%--I_tilde (pseudo-identity matrix used to seperate between streamfunction
%   points where pde is solved and those were Neumann BC is applied

    mats.I_til = get_I_til( parms ); 

%--

%--I_til_s2vort (used in Neumann BC)

    mats.I_til_s2vort = get_I_til_s2vort( parms ); 


%--

