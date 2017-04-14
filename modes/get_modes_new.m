clear all, close all, clc

%   Solves the generalized eigenvalue problem
%                       lambda ( Bx ) = Ax
%   using eigs,
%   where A is the Jacobian and B is the matrix [ RC 0; 0 0 ]

%NOTE: requires a base flow obtained from the file base.mat (you can use
%      the base_flow directory to compute it)

%--user specified info
    %# of evals and evects to compute
    k = 1;

    %specify where you want evals.
    lam_ref = 0 + 0.7 * i;
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


%--get initial guess

    addpath( '../base_flow/build_mats/')
%     v0 = get_v0( parms, mats, soln );

%--


%--Get A and B
    %Compute Jacobian (linearized about base flow)
    mats = assemble_Jac( parms, mats, soln );

    %Compute B
    mats = assemble_B( parms, mats, soln );
%--


%--compute the modes

    AmsigB = mats.A - lam_ref * mats.B;
    
%     condest( AmsigB ) 
%     
%     return

    [LL,UU,pp,qq,rr] = lu(AmsigB);

    %get matrix function handle for generalized e-val problem
    matfun = @(x) qq*(UU\(LL\(pp*(rr\ x ) ) ) );
    n = length(mats.A);
    
%     matfun = @(x) mats.A * x;

    opts.issym = 0;
    opts.isreal = 0;
    opts.disp = 2; %get full output from eigs
%     opts.p = 100;
    opts.tol = 1e-6;
    [V,D] = eigs( matfun, n, mats.B, k, lam_ref, opts);
    
    D = diag(D);
    imD = imag(D) ;
    imD = imD / (2*pi);
    D = real(D) + 1i * imD
    save( 'modes.mat', 'D','V', 'parms', 'mats' );

%--
