function K_flag = get_K_flag_FD4( Fint, xb, yb, chi, R_E, R_th, ...
            R_sh, beta0, h_el_0, parms, soln )
        
        
eps = 1e-6;
nb = length(xb);
nf = 2 * nb;

K_flag = zeros(3*nb);

for j = 1 : 3*nb
        
    %--get perturbation vars
        chi_p = zeros( 3*nb, 1);
        xbp = zeros(nf, 1);
        chi_p( j ) = 1;

    %--evaluate at x+h
        xbp( 1 : nb ) = parms.xb0( 1 : nb ) + soln.chi(1 : 3 : 3*nb-2)' + ...
                    eps * chi_p( 1 : 3 : 3*nb-2 )'; %x-body position
        xbp( nb + (1 : nb) ) = parms.xb0( nb + (1 : nb) ) + soln.chi(2 : 3 : 3*nb-1)' +... 
            eps * chi_p( 2 : 3 : 3*nb-1)';

        xb = xbp(1 : nb);
        yb = xbp( nb + (1 : nb ) );

        %compute updated geometry of flag relative to undeformed config.
        [c, s, h_el_d, qL, Fintph] = var_update(xb, yb, soln.chi + eps*chi_p, R_E, R_th, ...
        R_sh, beta0, h_el_0, parms);
    
    %--evaluate at x-h
        xbp( 1 : nb ) = parms.xb0( 1 : nb ) + soln.chi(1 : 3 : 3*nb-2)' - ...
                    eps * chi_p( 1 : 3 : 3*nb-2 )'; %x-body position
        xbp( nb + (1 : nb) ) = parms.xb0( nb + (1 : nb) ) + soln.chi(2 : 3 : 3*nb-1)' -... 
            eps * chi_p( 2 : 3 : 3*nb-1)';

        xb = xbp(1 : nb);
        yb = xbp( nb + (1 : nb ) );

        %compute updated geometry of flag relative to undeformed config.
        [c, s, h_el_d, qL, Fintmh] = var_update(xb, yb, soln.chi - eps*chi_p, R_E, R_th, ...
        R_sh, beta0, h_el_0, parms);
    
    %--evaluate at x+2h
        xbp( 1 : nb ) = parms.xb0( 1 : nb ) + soln.chi(1 : 3 : 3*nb-2)' + ...
                   2 * eps * chi_p( 1 : 3 : 3*nb-2 )'; %x-body position
        xbp( nb + (1 : nb) ) = parms.xb0( nb + (1 : nb) ) + soln.chi(2 : 3 : 3*nb-1)' +... 
           2 * eps * chi_p( 2 : 3 : 3*nb-1)';

        xb = xbp(1 : nb);
        yb = xbp( nb + (1 : nb ) );

        %compute updated geometry of flag relative to undeformed config.
        [c, s, h_el_d, qL, Fintp2h] = var_update(xb, yb, soln.chi + 2*eps*chi_p, R_E, R_th, ...
        R_sh, beta0, h_el_0, parms);
    
    %--evaluate at x-2h
        xbp( 1 : nb ) = parms.xb0( 1 : nb ) + soln.chi(1 : 3 : 3*nb-2)' - ...
                   2 * eps * chi_p( 1 : 3 : 3*nb-2 )'; %x-body position
        xbp( nb + (1 : nb) ) = parms.xb0( nb + (1 : nb) ) + soln.chi(2 : 3 : 3*nb-1)' -... 
           2 * eps * chi_p( 2 : 3 : 3*nb-1)';

        xb = xbp(1 : nb);
        yb = xbp( nb + (1 : nb ) );

        %compute updated geometry of flag relative to undeformed config.
        [c, s, h_el_d, qL, Fintm2h] = var_update(xb, yb, soln.chi - 2*eps*chi_p, R_E, R_th, ...
        R_sh, beta0, h_el_0, parms);
    
    
    
    %--add to column of stiffness matrix
        K_flag(:,j) = ( -Fintp2h + 8*Fintph - 8*Fintmh + Fintm2h) / (12*eps);
    
end  

%BCs

%clamped?
if parms.clamped == 'T'
   
    bc_type = [1 1 1];
    
else
    
    bc_type = [1 1 0];
    
end

count = 0;

if parms.inverted == 'T'
   
    for j = 3*nb - 2 : 3*nb
       
        count = count + 1;
        
        if bc_type(count) == 1
            K_flag(j,:) = 0;
%             K_flag(:,j) = 0;
            K_flag(j,j) = 1;
        end
        
    end
    
else
    
    for j = 1 : 3
       
        count = count + 1;
        
        if bc_type(count) == 1
            K_flag(j,:) = 0;
%             K_flag(:,j) = 0;
            K_flag(j,j) = 1;
        end
        
    end
    
    
end


% % % cond( K_flag )
% % % 
% % % pause

