clear all, close all, clc


%%
addpath('../build_mats/')

% Contour maximum values and number of contour levels
% Vorticity
cmax_w = 3;
clev_w = 20;
clevs = linspace( -cmax_w, cmax_w, clev_w );

% Range for plots
range = [-1 5 -2 2];

load('cmap.mat')

for it = [2000 : 100 : 2000]

    load(['outputs/runvars_it_',num2str(it),'.mat'])

    gamma = soln.gamma;
    
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

            omega = gamma( 1 : (parms.m-1)*(parms.n-1), lev ) / delta^2;

            omega(omega > cmax_w ) = cmax_w;
            omega(omega < -cmax_w ) = -cmax_w;

            Omega(:,:,lev) = transpose( reshape( omega, parms.m-1, parms.n-1 ) );                                  

   
        %--


    end
    
    %plot them
    
    figure(1), clf

    for j = parms.mg : -1 : 1
        
        figure(1), hold on
        contourf(Xv(:,:,j), Yv(:,:,j), Omega(:,:,j), clevs, ...
            'edgecolor','none'); shading flat;

        colormap( cmap )
        axis equal
%         axis(range)
             
    end
    
    %plot body
    fill(parms.xb( 1 : parms.nb ), parms.xb( 1+parms.nb : 2*parms.nb ),'k'  )

    pause
    
    
end
