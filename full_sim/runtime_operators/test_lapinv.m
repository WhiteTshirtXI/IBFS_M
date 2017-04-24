clear all, close all, clc

addpath('./build_mats/')

parms.m = 10; parms.n = 10; parms.mg = 1;

m = parms.m; n = parms.n;


ngam = (parms.m-1) * (parms.n-1);

display('getting curl and rot matrices')    
[mats.C, mats.R] = get_curls( parms ); 

mats.RC = mats.R * mats.C;

% display('getting LU factorization of RC. Patience, grasshoppa...')
% [mats.LRC,mats.URC,mats.pRC,mats.qRC,mats.rRC] = lu(mats.RC);
% 
% %get matrix function handle for solving linear system
% mats.invRC = @(x) mats.qRC*(mats.URC\(mats.LRC\(mats.pRC*(mats.rRC \ x ) ) ) );

b = randn( ngam, 1) ;

Lap = -mats.RC; 

tic
% x1 = -mats.invRC( b );
x1 = Lap \ b;
toc


ii = 1 : (m-1); jj = 1 : (n-1);
[ii, jj] = meshgrid( ii, jj );
mats.lam = 2 .* ( cos( pi .* ii ./ m ) + cos( pi .* jj ./ n ) - 2);

tic
x2 = lapinv( b,  parms, mats );
toc


max(max( abs( x1 - x2 ) ) )






