function [] = main_fun( soln_lin, soln, mats, parms )

addpath('./build_mats/')

   
%---advance in time

    for it = parms.it_start : parms.it_stop
        
        %output occasionally to tell us we're advancing in time
        if mod( it, 10 ) == 0 & it > 0
            display( ['Advancing to time step ', num2str( it+1 )] )
        end
        
       
        %advance a time step and return circulation (gamma), vel flux (q), and 
        %surface stress (fb)
        [soln_lin, mats] = advance_lin( it, soln_lin, parms, mats, soln );

        %save other variables if at a save point
        if ( ( mod( it, parms.it_save ) == 0 ) | ( it == parms.it_stop ) )

            save(['outputs/runvars_it_',num2str(it),'.mat'], ...
                'soln_lin','soln','mats', 'parms' );
        end
        
    end
    
%---





