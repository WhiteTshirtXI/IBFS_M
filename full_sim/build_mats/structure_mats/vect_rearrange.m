function vect = vect_rearrange( vect )

%Rearrange element matrices to correspond to the following arrangement of
%unknowns and eqns:
%           [u1, w1, theta1, u2, w2, theta2]
%(the matrices are originally arranged as [u1, u2, w1, theta1, w2, theta2]

v_store = vect;

vect(2) = v_store(3);
vect(3) = v_store(4);
vect(4) = v_store(2);