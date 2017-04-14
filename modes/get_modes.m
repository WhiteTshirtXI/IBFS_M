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
    lam_ref = 0.1 + 0.37 * 2*pi * 1i;
%--

%--check to make sure we have a base flow

    file_base = '../base_flow/base.mat';
%     file_base = 'base_rigid.mat';

    %don't have a base flow
    if exist( file_base , 'file') ~= 2
        error( ['Base flow required! Run get_base.m in the '...
            ' base_flow directory to compute it'])
    else
        load(file_base)
        %Parameters for flag
        parms.R_rho = 30;
        parms.R_E = 1e3;
        parms.R_sh = 1e-7;
        parms.R_th = 0.01;
        parms.inverted = 'F';
        parms.clamped = 'F';
% %         parms.xb0 = soln.xb;
    end

%--

addpath( '../base_flow/')
addpath( '../base_flow/build_mats/')

    %first build the constituent matrices
    mats = get_mats_preproc( parms );

    mats = get_mats( parms, mats, soln );

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
    opts.tol = 1e-10;
    [V,D] = eigs( matfun, n, mats.B, k, lam_ref, opts);
    
    D = diag(D);
    imD = imag(D) ;
    imD = imD / (2*pi);
    D = real(D) + 1i * imD
    save( 'modes.mat', 'D','V', 'parms', 'mats' );

%--
