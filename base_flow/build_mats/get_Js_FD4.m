function [J_Ef, J_Wf, J_Es] = get_Js_FD4( parms, mats, soln )

%build matrices corresponding to third rank tensor contributions from
%matrices that depend on body position (E, ET, W)

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

eps = 1e-6; %A small number for the finite difference approximation.


fb_base = soln.fb;
s_base = soln.s;


%Get reference E, ET, W at current position
ET_ref = mats.ET;
E_ref = mats.E ;
W_flag_ref = mats.W_flag;


%build the Js 1 column at a time
J_Ef = zeros( nq, nf );
J_Wf = zeros( nf, nf );
J_Es = zeros( nf, nq );

for j = 1 : nf

    %perturb body position
    xbp = zeros( 1, 2* nb );
    xbp(j) = 1;

    %--evaluate at x+h
       
        %Get E, ET, W at perturbed position
        ET_ph = get_ET( soln.xb + eps * xbp , parms, soln );
        E_ph = ET_ph' ;
        W_flag_ph = get_W_flag( ET_ph, E_ph, parms );
        
        
    %--
        
    %--evaluate at x-h
    
        %Get E, ET, W at perturbed position
        ET_mh = get_ET( soln.xb - eps * xbp , parms, soln );
        E_mh = ET_mh' ;
        W_flag_mh = get_W_flag( ET_mh, E_mh, parms );
    
    %--
    
    %--evaluate at x+2h
    
        %Get E, ET, W at perturbed position
        ET_p2h = get_ET( soln.xb + 2 * eps * xbp , parms, soln );
        E_p2h = ET_p2h' ;
        W_flag_p2h = get_W_flag( ET_p2h, E_p2h, parms );
    
    %--
    
    %--evaluate at x-2h
    
        %Get E, ET, W at perturbed position
        ET_m2h = get_ET( soln.xb - 2 * eps * xbp , parms, soln );
        E_m2h = ET_m2h' ;
        W_flag_m2h = get_W_flag( ET_m2h, E_m2h, parms );
    
    %--
    
    %--
    

        J_Ef(:,j) = ( -ET_p2h * fb_base + 8*ET_ph * fb_base - ...
            8*ET_mh*fb_base + ET_m2h*fb_base) / (12*eps);
        
        
        
        J_Wf(:,j) = ( -W_flag_p2h * fb_base + 8*W_flag_ph * fb_base -...
            8*W_flag_mh * fb_base + W_flag_m2h* fb_base ) / (12*eps);
        
        
        
        J_Es(:,j) = ( -E_p2h * mats.M_vel * ( mats.C * s_base + parms.q0 ) + ...
           8*E_ph * mats.M_vel * ( mats.C * s_base + parms.q0 ) -...
           8*E_mh* mats.M_vel * ( mats.C * s_base + parms.q0 ) + ...
           E_m2h * mats.M_vel * ( mats.C * s_base + parms.q0 ) ) / (12*eps);
    
end


J_Ef = mats.R * J_Ef;
J_Wf = mats.Q_flag * J_Wf * mats.Itilde_flag;





