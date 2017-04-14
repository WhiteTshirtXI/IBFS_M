function [] = main_fun( parms, mats, soln )

addpath('./build_mats/')

%--Double-check that m and n were specified correctly
    if ( mod(parms.m, 4) ~=0 | mod(parms.n, 4) ~= 0 )
        error( 'Error: parms.m and parms.n must be divisible by 4')
    end

%--

%---Build undeformed configuration for body
    
    h = parms.len / parms.m;
    parms.ds = 2*h; %IB grid spacing is 2 times the fluid grid spacing
    
    xbv = 0 : parms.ds : 1;
    ybv = zeros(size(xbv));
    parms.nb = length(xbv);
  
    parms.xb0 = zeros( 1, 2*parms.nb);
    parms.xb0( 1 : parms.nb ) = xbv;
    parms.xb0( parms.nb + 1 : 2*parms.nb ) = ybv;
    
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
    nb = parms.nb;
    nf = 2*parms.nb;
    %grid spacing on finest grid (used a lot)
    h = parms.len / parms.m;
%--
    
%---preprocessing: 

    %build matrices if they haven't been build
%     if exist( file_start , 'file') ~= 2
        tic
        %build and store matrices using sparse operations
        display('------------------------------------------------------------')
        display('Pre-processing stage: building and storing matrices for run')
        
        %first build the constituent matrices
        mats = get_mats_preproc( parms );
        
        pre_time = toc;

        display('Done building and storing matrices')
        display(['Spent ',num2str(pre_time),' secs on pre-processing.']) 
        display('------------------------------------------------------------')
%     end
    
    file_start = 'base.mat';

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
        soln.s = zeros( ngam, 1); %Streamfcn
        soln.fb = zeros( nf, 1); %Surface stress
        soln.chi = zeros( 3*nb, 1); %IB position (includes rotations)
% % %         soln.chi( 2 : 3 : 3*nb-1) = 0.1 * linspace(1,0, nb);
        soln.xb =  parms.xb0; %also store positions w/out rotations
% % %         soln.xb( nb + (1:nb) ) = soln.xb( nb + (1:nb) ) + soln.chi( 2 : 3 : 3*nb -1 )';
% % %         xbshow = soln.xb
        soln.zeta = soln.chi; %IB velocity
    end
    
% % %      soln.chi = zeros( 3*nb, 1); %IB position (includes rotations)
% % %      soln.xb =  parms.xb0; %also store positions w/out rotations
% % %      soln.zeta = soln.chi; %IB velocity

    y = [soln.s; soln.zeta; soln.chi; soln.fb];

    iter = 1;
    
    while err > parms.tol
        
        display('      ------------------------------------')
        display(['      Beginning iteration ', num2str( iter )])
        
        
        %compute matrices that change with change in state
        mats = get_mats( parms, mats, soln );
        
        
        %Update Jacobian
        mats = assemble_mats( parms, mats, soln );
        
        
        [LL,UU,pp,qq,rr] = lu(mats.Dr);

        
        %Update function eval
        rold = r_eval( parms, mats, soln );
        
        if iter == 1
            
            err = norm( rold );
        else
            
            err = norm( rold ) / norm( y );
        end
        display(['      Error = ', num2str( err )])
        
        if err > parms.tol
            
            display(['      Updating guess... ', num2str( iter )])
            display('      Solving linear system for dy... ')
            dy = -qq*(UU\(LL\(pp*(rr\ rold ) ) ) );
            
            
            dyinf = max(abs( dy ) )
            
            
            y = y + dy;
            soln.s = y(1 : ngam);
            soln.zeta = y( ngam + (1 : 3*nb ) ); 
            soln.chi = y( ngam + 3*nb + (1 : 3*nb) );
            soln.fb = y( ngam + 6*nb + (1 : nf ) );
            
                        
            soln.q = mats.C * soln.s; %get vel flux
            soln.gamma = mats.R * soln.q; %get circulation
            
           
            soln.xb( 1 : nb ) = parms.xb0( 1 : nb ) + ...
                soln.chi( 1 : 3 : 3*nb-2 )'; %x-body position
            soln.xb( nb + (1 : nb) ) = parms.xb0( nb + (1 : nb) ) + ... 
                soln.chi( 2 : 3 : 3*nb-1)'; %y-body position
                    
        end
        
        display('      ------------------------------------')

        iter = iter + 1;
        save('base.mat', 'soln','mats', 'parms' );
        
    end
    
    display('------------------------------------------------------------')

%---





