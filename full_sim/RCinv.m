function x = RCinv( b, parms, mats )

%Solve RC * x = b, where -RC is the standard 5-point stencil Laplacian.

m = parms.m; n = parms.n;

%reshape for inversion in fourier space
b = reshape( b, parms.m-1, parms.n-1 ); 

%Solve by transforming to and from Fourier space and scaling by evals   
scl = 4 / (m * n); %scale that makes the inverse and direct transforms equal

lam_sc = -mats.lam ./ scl; %The -sign is because RC = -L

x = dst( transpose( dst( dst( transpose( dst( b ) ) ) ./ lam_sc ) ) );

%give output in same size as input b (before being reshaped)
ngam = (parms.m-1) * (parms.n-1);
x = reshape( x, ngam, 1);
    