function x = lapinv( b, parms, mats )

%reshape for inversion in fourier space
b = reshape( b, parms.m-1, parms.n-1); 

%Solve by transforming to and from Fourier space and scaling by evals   
x = idst( transpose( idst( dst( transpose( dst( b ) ) ) ./ mats.lam ) ) );

%give output in same size as input b (before being reshaped)
ngam = (parms.m-1) * (parms.n-1);
x = reshape( x, ngam, 1);
    