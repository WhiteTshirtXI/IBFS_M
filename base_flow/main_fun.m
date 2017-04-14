function [] = main_fun( parms, mats, soln )

addpath('./build_mats/')

%--Double-check that m and n were specified correctly
    if ( mod(parms.m, 4) ~=0 | mod(parms.n, 4) ~= 0 )
        error( 'Error: parms.m and parms.n must be divisible by 4')
    end

%--

%---Load body

    if parms.body == 'cyl'
        
        %radius of body is 1/2
        [xb, ds] = build_cylinder( 1/2, parms.len/parms.m );
    else
        Warning(['Requested body type is not supported. \n'...
            'Must provide your own body points as a vector xb.'])
    end
    
    load( 'body.mat' );
    parms.xb = xb;
    parms.nb = length( xb ) / 2;
    parms.ds = ds;
        
%---


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
    nf = 2*parms.nb;
    %grid spacing on finest grid (used a lot)
    h = parms.len / parms.m;
%--
    
%---preprocessing: 

    %build matrices if they haven't been build
    file_start = 'base.mat';
%     if exist( file_start , 'file') ~= 2
        tic
        %build and store matrices using sparse operations
        display('------------------------------------------------------------')
        display('Pre-processing stage: building and storing matrices for run')
        
        %first build the constituent matrices
        mats = get_mats( parms );
        
        pre_time = toc;

        display('Done building and storing matrices')
        display(['Spent ',num2str(pre_time),' secs on pre-processing.']) 
        display('------------------------------------------------------------')
%     end
    
    %build and store q0 if not loading file
    if exist( file_start , 'file') ~= 2
        
        q0 = zeros( nq, 1 );
        for j = 1 : parms.mg

            %grid spacing on current grid
            hc = h * 2^(j-1);

            %index of x-vel flux points for current grid
            if j == 1
                ind_x = 1 : get_velx_ind( parms.m-1, parms.n, j, parms );
            else
                ind_x = 1 + get_velx_ind( parms.m-1, parms.n, j-1, parms ) : ...
                    get_velx_ind( parms.m-1, parms.n, j, parms );
            end

            %write fluid velocity flux in body-fixed frame
            q0( ind_x ) = -parms.U_body * hc;
        end
        parms.q0 = q0;
        
    end
%---


%---N-R loop

    display('------------------------------------------------------------')
    display('Entering Newton-Raphson loop' )

    err = 1;
    
    %create initial guess if not loaded from file
    if exist( file_start , 'file') ~= 2
        soln.s = zeros( ngam, 1);
        soln.fb = zeros( nf, 1);
    end
    

    y = [soln.s; soln.fb];

    iter = 1;
    
    while err > parms.tol
        
        display('      ------------------------------------')
        display(['      Beginning iteration ', num2str( iter )])
        
        %Update Jacobian
        mats = assemble_mats( parms, mats, soln );
        
        
        [LL,UU,pp,qq,rr] = lu(mats.Dg);

        
        %Update function eval
        gold = g_eval( parms, mats, soln );
        
        if iter == 1
            
            err = norm( gold );
        else
            goldsm = gold( 1 : get_vort_ind( parms.m-1, ...
                parms.n-1, 2, parms ) );
            ysm = y( 1 : get_vort_ind( parms.m-1, ...
                parms.n-1, 2, parms ) );

            err = norm( goldsm ) / norm( ysm );
        end
        display(['      Error = ', num2str( err )])
        
        if err > parms.tol
            
            display(['      Updating guess... ', num2str( iter )])
            display('      Solving linear system for dy... ')
            dy = -qq*(UU\(LL\(pp*(rr\ gold ) ) ) );
            
            
%             dyinf = max(abs( dy ) )
            
            
            y = y + dy;
            soln.s = y(1 : ngam);
            soln.fb = y( ngam + 1 : end );
                        
            soln.q = mats.C * soln.s; %get vel flux
            soln.gamma = mats.R * soln.q; %get circulation
        
        end
        
        display('      ------------------------------------')

        iter = iter + 1;
        save('base.mat', 'soln','mats', 'parms' );
        
    end
    
    display('------------------------------------------------------------')

%---





