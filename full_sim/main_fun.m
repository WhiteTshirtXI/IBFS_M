function [] = main_fun( soln, mats, parms )

addpath('./build_mats/')

%--Double-check that m and n were specified correctly
    if ( mod(parms.m, 4) ~=0 | mod(parms.n, 4) ~= 0 )
        error( 'Error: parms.m and parms.n must be divisible by 4')
    end

%--

%---Build flag

    h = parms.len / parms.m;
    parms.ds = 2*h;
    
    %Reference configuration for flag
    xb0x = 0 : parms.ds : 1;
    xb0y = zeros( size( xb0x ) );
    
    parms.xb0 = [xb0x xb0y];
    parms.nb = length( xb0x );
    soln.xb = parms.xb0; %body position begins in undeformed state
    
%---
    
%---preprocessing: 

    if isempty( mats )
        tic
        %build and store matrices using sparse operations
        display('------------------------------------------------------------')
        display('Pre-processing stage: building and storing matrices for run')
        mats = get_mats_preproc( parms );

        pre_time = toc;

        display('Done building and storing matrices')
        display(['Spent ',num2str(pre_time),' secs on pre-processing.']) 
        display('------------------------------------------------------------')
    end
    
%---
    

%---advance in time

    for it = parms.it_start : parms.it_stop
        
        %output occasionally to tell us we're advancing in time
        if mod( it, 1 ) == 0 & it > 0
            display( ['Advancing to time step ', num2str( it+1 )] )
        end
        
        %advance a time step and return circulation (gamma), vel flux (q), and 
        %surface stress (fb)
        [soln,parms,mats] = advance( it, parms, mats, soln );
        
        
        %save other variables if at a save point
        if ( ( mod( it, parms.it_save ) == 0 ) | ( it == parms.it_stop ) )

            save(['outputs/runvars_it_',num2str(it),'.mat'], ...
                'soln', 'parms' );
        end
        
    end
    
%---





