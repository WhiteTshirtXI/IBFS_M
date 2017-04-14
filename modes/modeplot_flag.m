clear all, close all, clc


%%
addpath('../base_flow/build_mats/')

% Contour maximum values and number of contour levels
% Vorticity
% cmax_w = 3;
% clev_w = 20;
% clevs = linspace( -cmax_w, cmax_w, clev_w );

% Range for plots
range = [-5 10 -5 5];

load('cmap.mat')

load('../rigid_base/base.mat')

load('modes.mat')

k = 1;

%--Various variables 
    %# of x-vel (flux) points
    nu = get_velx_ind( parms.m-1, parms.n, parms.mg, parms );
    %# of y-vel (flux) points
    nv = get_vely_ind( parms.m, parms.n-1, parms.mg, parms );
    %Total # of vel (flux) points
    nq = nu + nv;
    %# of vort (circ) points
    ngam = get_vort_ind( parms.m-1, parms.n-1, parms.mg, parms );
    %# of body points
    nb = parms.nb;
    %# of surface stress points
    nf = nb * 2;
%--

for jj = 1 : 2

    if jj == 1
        Vv = real( V( ngam + 3*nb + ( 1 : 3*nb ), k) );
        str = 'Re($\chi_p$)';
    else
        Vv = imag( V( ngam + 3*nb + ( 1 : 3*nb ), k) );
        str = 'Im($\chi_p$)';
    end
        
    Vvy = Vv(2 : 3 : 3*nb-1)';
    Vvx = Vv( 1 : 3 : 3*nb-2)';
    
    xbx = parms.xb0( 1 : nb ); xby = parms.xb0( nb + (1 : nb ) );
    
    xbp = xbx + Vvx; ybp = xby + Vvy;


    figure(10), subplot( 2,1, jj)
    plot( xbp, ybp, 'k','linewidth',2 )
%     axis equal
    axis([0 1 -0.005 0.005])
    set(gca,'TickLabelInterpreter','latex','fontsize',14)
    title( str, 'fontsize', 18, 'interpreter','latex')
end

