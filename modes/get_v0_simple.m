function v0 = get_v0_simple( parms, mats, soln )

%initial guess is random noise advanced with the linear solver ten steps in
%time

if exist( 'v0.mat', 'file' ) ~= 2 
    
        parms.dt = 1e-3;

    %--define number of points for various quantities

        %# of x-vel (flux) points
        nu = get_velx_ind( parms.m-1, parms.n, parms.mg, parms );
        %# of y-vel (flux) points
        nv = get_vely_ind( parms.m, parms.n-1, parms.mg, parms );
        %Total # of vel (flux) points
        nq = nu + nv;
        %# of vort (circ) points
        ngam = get_vort_ind( parms.m, parms.n-1, parms.mg, parms );
        %# of streamfunction points (same as # of vort points)
        nstrm = ngam;
        %# of body points
        nb = parms.nb;
        %# of surface stress points
        nf = nb * 2;
        %time step size in physical time
        dt = parms.dt;
        %Grid spacing on IB
        ds = parms.ds;

    %--

    % % % %--generate random noise on all but largest grid
        %# of vort (circ) points on all but largest grid
        ngam_mgm1 = get_vort_ind( parms.m-1, parms.n-1, parms.mg-1, parms );
        %# of surface stress points
        nf = parms.nb * 2;
        v0 = zeros( ngam + nf, 1 );
        v0( 1 : ngam_mgm1 ) = randn( ngam_mgm1, 1 );
        v0( ngam + 1 : ngam + nf ) = randn( nf, 1 );
        
    % % % %--
    
    save( 'v0.mat', 'v0' )

else
    
    load( 'v0.mat')
end
