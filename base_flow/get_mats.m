function mats = get_mats( parms, mats, soln )

%Build and store matrices required for simulation


%!!!!!!!Operators for fluid
    
    %---Smearing and interpolation matrices

        display('getting smearing and interpolation matrices')
        parms.supp = 6; %currently using a delta function with 
                        %a support of 6 grid cells...


        mats.ET = get_ET( soln.xb, parms, soln );

        mats.E = mats.ET' ;

    %---
    
%!!!!!!


%!!!!!!Operators for flag

    %define variables for ease of use
    nb = parms.nb; 
    xb0 = parms.xb0(1:nb); 
    yb0 = parms.xb0(nb + (1:nb));
    xb = soln.xb(1:nb); 
    yb = soln.xb(nb + (1:nb));
    chi = soln.chi;
    R_E = parms.R_E;
    R_th = parms.R_th;
    R_sh = parms.R_sh;
    R_rho = parms.R_rho;
    nel = nb - 1;
    
    %---mass matrix for flag
    
        %compute flag geometry in reference geometry (usually undeformed)
        [beta0, h_el_0] = var_init( xb0, yb0 );
        
        h_el_0 = h_el_0(1);
    
        mats.M_flag = get_M_flag( R_rho, R_th, h_el_0, nel, nb, parms );
    
    %---
    
    %---Q matrix for flag
    
        mats.Q_flag = get_Q_flag( nb, h_el_0, parms );
    
    %---
    
    %---Get W matrix for flag
    
        mats.W_flag = get_W_flag( mats.ET, mats.E, parms );
        
    %---
    
    %---Itilde matrix for flag (takes a 3*nb vector and converts it to a
    %   2*nb vector)
    
        mats.Itilde_flag = get_Itilde_flag( nb );
    
    %---
    
    
    %---Stiffness matrix for flag
    
        %compute updated geometry of flag relative to undeformed config.
        [c, s, h_el_d, qL, mats.Fint] = var_update(xb, yb, chi, R_E, R_th, ...
            R_sh, beta0, h_el_0, parms);

%         mats.Fint
        
    
        %Update stiffness matrix given this new body position
        mats.K_flag_rot = get_K_flag( R_E, R_th, R_sh, c, s, h_el_d, parms ) ;
%         mats.K_flag = get_K_flag_rig( R_E, R_th, R_sh, h_el_0, ...
%             c, s, h_el_d, qL, parms ) ;
% 
          mats.K_flag = get_K_flag_FD( mats.Fint, xb, yb, chi, R_E, R_th, ...
            R_sh, beta0, h_el_0, parms, soln );

    %---                                   
    
    %---Contributions from third rank tensors
    
        [mats.J_Ef, mats.J_Wf, mats.J_Es] = get_Js( parms, mats, soln );
    
%         [mats.J_Ef, mats.J_Wf, mats.J_Es] = get_Js_FD4( parms, mats, soln );

        
    %---

   
%!!!!!!



