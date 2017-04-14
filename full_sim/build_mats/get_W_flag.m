function W_flag = get_W_flag( ET, E, parms )

%W_flag = E * W * ET, where W is a weighting matrix that preserves the
%integral weight of the surface stress.

nf = length( ET( 1,:) ) ;
nq = length( ET( :,1) ) ;

%First build weighting matrix W:

%vector of ones
one_v = ones( nf, 1);

wght = ET * one_v;

wght( abs(wght) <= 1e-10 ) = 0;

wght( wght ~= 0 ) = 1 ./ wght( wght ~= 0 );

Wght = sparse( 1 : nq, 1 : nq, wght ) ;

W_flag = E * Wght * ET;


% %scaling to account for different grid spacing
% h = parms.len / parms.m;
% 
% W_flag = W_flag * (h / parms.ds);
