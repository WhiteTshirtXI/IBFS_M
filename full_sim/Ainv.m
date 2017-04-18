function x = Ainv( b, lev, parms, mats )

%reshape for inversion in fourier space
b = reshape( b, parms.m-1, parms.n-1); 

%Solve by transforming to and from Fourier space and scaling by evals 
hc = (parms.len / parms.m) * 2^( lev - 1);
lamtil = 1 + mats.lam * parms.dt/( 2 * parms.Re * hc^2 );
x = idst( transpose( idst( dst( transpose( dst( b ) ) ) ./ lamtil ) ) );

%give output in same size as input b (before being reshaped)
ngam = (parms.m-1) * (parms.n-1);
x = reshape( x, ngam, 1);
    