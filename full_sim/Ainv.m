function x = Ainv( b, lev, parms, mats )

%Solve (I + dt/2 * Beta * RC) * x = b for x, where Beta = 1/(Re * h^2) 

%reshape for inversion in fourier space
b = reshape( b, parms.m-1, parms.n-1); 

%Solve by transforming to and from Fourier space and scaling by evals 
hc = (parms.len / parms.m) * 2^( lev - 1);
lamtil = 1 + mats.lam * parms.dt/( 2 * parms.Re * hc^2 );

%scale to make inverse and direct transforms equal
scl = 4 / (parms.m * parms.n ); 
x = dst( transpose( dst( dst( transpose( dst( b ) ) ) ./ lamtil ) ) ) .* scl;

%give output in same size as input b (before being reshaped)
ngam = (parms.m-1) * (parms.n-1);
x = reshape( x, ngam, 1);
    