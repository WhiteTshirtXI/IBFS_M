clear all, close all, clc


%%
addpath('./build_mats/')

% Contour maximum values and number of contour levels
% Vorticity
cmax_w = 3;
clev_w = 20;
clevs = linspace( -cmax_w, cmax_w, clev_w );

% Range for plots
range = [-5 10 -5 5];

load('cmap.mat')


    load(['base.mat'])

    gamma = soln.gamma;
    
    s = (mats.RC ) \ gamma;
    
    Xv = zeros( parms.n-1, parms.m-1, parms.mg );
    Yv = Xv;
    Omega = Xv;

    %Get x-y points and vorticity
    for lev = 1 : parms.mg

        fac = 2^(lev-1);

        % Grid spacing in both directions for current grid
        delta = parms.len ./ parms.m *fac;

        % Offset in x direction for current grid
        offx = 2^(lev-1) * parms.len/2 - parms.len/2 + parms.offx ;

        % Offset in y direction for current grid
        offy = 2^(lev-1) * (parms.n*parms.len/parms.m)/2 - ...
        (parms.n*parms.len/parms.m)/2 + parms.offy ;


        %--get grid points


            xv = delta : delta : (parms.m-1) * delta;
            yv = delta : delta : (parms.n-1) * delta;

            xv = xv - offx;
            yv = yv - offy;

            [Xv(:,:,lev), Yv(:,:,lev)] = meshgrid( xv, yv );
 
        %--
    
        %--get vorticity

        if lev == 1

            ind_s = 1;
            ind_e = get_vort_ind( parms.m-1, parms.n-1, lev, parms );
            omega = gamma(ind_s : ind_e ) / delta^2;

            omega(omega > cmax_w ) = cmax_w;
            omega(omega < -cmax_w ) = -cmax_w;

            Omega(:,:,lev) = transpose( reshape( omega, parms.m-1, parms.n-1 ) );                                  

            %store omega for coarse-grid interpolation
            omega_s = omega;
        else
            
            omega = zeros( (parms.m-1) * (parms.n-1), 1 );
            
            %--need some fancy indexing for what follows:
            
                %get indices of current grid (indices for overlapping
                %region are zero)
                indvortx = repmat(1 : parms.m-1,[1,parms.n-1]);
                indvorty = repelem(1 : parms.n-1, parms.m-1);
                vort_ind = get_vort_ind(indvortx,indvorty,lev,parms);
            
                %get all indices on gridlevel 1
                indvortx = repmat(1 : parms.m-1,[1,parms.n-1]);
                indvorty = repelem(1 : parms.n-1, parms.m-1);
                vort_ind_f = get_vort_ind(indvortx,indvorty,1,parms);
                
                %Get indices on gridlevel 1 just for overlapping part
                indvortx = repmat(2 : 2 : parms.m-2,[1,parms.n/2-1]);
                indvorty = repelem(2 : 2 : parms.n-2, parms.m/2-1);
                vort_ind_s = get_vort_ind(indvortx,indvorty, 1, parms);
                
                %indexing on the vector gamma starts at 1...
                n_sub = get_vort_ind( parms.m-1, parms.n-1, lev-1, parms);
                
                %non-overlapping indices:
                ind_no_over = ( (vort_ind - n_sub) > 0  );
                ind_no_over = vort_ind_f( ind_no_over ~= 0 );
                
                %overlapping indices:
                ind_over = ( (vort_ind - n_sub) < 0  ); 
                ind_over = vort_ind_f( ind_over ~= 0 );
            
            %--
            
            
            
            %Get part that doesn't overlap with finer grid
            omega( ind_no_over ) = gamma( vort_ind( vort_ind ~= 0 ) )/ delta^2;
            
            
%             %for overlapping part...
            omega( ind_over ) = omega_s( vort_ind_s );
%                         
            Omega(:,:,lev) = transpose( reshape( omega, parms.m-1, parms.n-1 ) );
            
            %save omega for next coarse grid
            omega_s = omega;

        end

        %--



    end
    
    %plot them
    
    figure(1), clf

    for j = parms.mg : -1 : 1
        
        figure(1), hold on
%         contourf(Xv(:,:,j), Yv(:,:,j), Omega(:,:,j), clevs, ...
%             'edgecolor','none'); shading flat;

        pcolor(Xv(:,:,j), Yv(:,:,j), Omega(:,:,j) ); shading interp;

        colormap( cmap )
        axis equal
%         axis(range)
                
    end
    
    %plot body
    fill(parms.xb( 1 : parms.nb ), parms.xb( 1+parms.nb : 2*parms.nb ),'k'  )


