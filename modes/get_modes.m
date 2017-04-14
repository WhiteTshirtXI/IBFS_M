clear all, close all, clc

%   Solves the generalized eigenvalue problem
%                       lambda ( Bx ) = Ax
%   using eigs,
%   where A is the Jacobian and B is the matrix [ RC 0; 0 0 ]

%NOTE: requires a base flow obtained from the file base.mat (you can use
%      the base_flow directory to compute it)

%--user specified info
    %# of evals and evects to compute
    k = 6;

    %specify where you want evals (0 is a good default).
    lam_ref = 0;
%--

%--get initial guess

% % %     v0 = get_v0( parms, mats, soln );

    load('../lin_sim/outputs/runvars_it_26000.mat')
    
    q = soln_lin.q;  gam = soln_lin.gamma;
    
    
    %# of vort (circ) points
    ngam = get_vort_ind( parms.m-1, parms.n-1, parms.mg, parms );
   
    %# of surface stress points
    nf = parms.nb * 2;
    
    v0 = zeros( ngam + nf, 1);
    v0( 1 : ngam ) = ( mats.RC ) \ gam;
    v0( ngam + 1 : ngam + nf ) = soln_lin.fb * parms.ds / (parms.len/parms.m );
     
    opts.v0 = v0;
    
    clear soln parms mats soln_lin
%--


%--check to make sure we have a base flow

    file_base = '../base_flow/base.mat';

    %don't have a base flow
    if exist( file_base , 'file') ~= 2
        error( ['Base flow required! Run get_base.m in the '...
            ' base_flow directory to compute it'])
    else
        load(file_base)
    end

%--


%--Get A and B
    addpath( '../base_flow/build_mats/')
    %Compute Jacobian (linearized about base flow)
    mats = assemble_Jac( parms, mats, soln );

    %Compute B
    mats = assemble_B( parms, mats, soln );
%--


%--compute the modes

            c = 1 / 0.7i;
            n = length( mats.A );
%             x = ones( n, 1 );
            Afun = @(x) ( mats.B - c * mats.A ) \ ( mats.A * x );

% % %     %get matrix function handle for generalized e-val problem
% % %     Afun = @(x) mats.A \ ( mats.B * x );
% % %     n = length(mats.A);
% % %     c = 0;
    
    opts.issym = 0;
    opts.isreal = 0;
    opts.disp = 2; %get full output from eigs
%     opts.p = 100;
    opts.tol = 1e-6;
    [V,D] = eigs( Afun, n, k, c, opts);
    
    D = 1./diag(D);
    imD = imag(D) ;
    imD = imD / (2*pi);
    D = real(D) + 1i * imD
    save( 'modes.mat', 'D','V', 'parms', 'mats' );

%--
