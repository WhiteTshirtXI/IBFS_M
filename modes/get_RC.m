function mats = get_RC( parms )

%Build and store matrices required for simulation


    %---Curl matrices

        display('getting curl and rot matrices')    
        [mats.C, mats.R] = get_curls( parms ); 

        mats.RC = mats.R * mats.C;

    %---

 