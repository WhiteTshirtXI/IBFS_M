clear all, close all, clc


%%
addpath('../build_mats/')

% Contour maximum values and number of contour levels
% Vorticity
cmax_v = 2;
clev_v = 20;
clevs = linspace( 0, cmax_v, clev_v );

% Range for plots
range = [-4 10 -5 5];

load('cmap.mat')

for it = 1000 : 1000 : 10000

    load(['outputs/runvars_it_',num2str(it),'.mat'])

    q = soln.q + soln.q0;

%     q = soln.q;
    
    Xv = zeros( parms.n, parms.m-1, parms.mg );
    Yv = Xv;
    U = Xv;

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
            yv = delta/2 : delta : (parms.n-1/2) * delta;

            xv = xv - offx;
            yv = yv - offy;

            [Xv(:,:,lev), Yv(:,:,lev)] = meshgrid( xv, yv );
 
        %--
    
        %--get velocity

            velu = q( 1 : (parms.m-1)*parms.n, lev ) / delta;

%             velu(velu > cmax_v ) = cmax_v;
%             velu(velu < 0 ) = 0;

            U(:,:,lev) = transpose( reshape( velu, parms.m-1, parms.n ) );                                  

        %--



    end
    
    %plot them
    
    figure(1), clf

    for j = 2 : -1 : 1
        
        figure(1), hold on
        contourf(Xv(:,:,j), Yv(:,:,j), U(:,:,j), clevs, ...
            'edgecolor','none'); shading flat;

        colormap( cmap )
        axis equal
%         axis(range)
                
    end
    
    %plot body
    fill(parms.xb( 1 : parms.nb ), parms.xb( 1+parms.nb : 2*parms.nb ),'k'  )


    pause
end
