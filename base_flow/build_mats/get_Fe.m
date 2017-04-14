function F_e = get_Fe( f, q, h )

%Get element-wise rhs vector F_e (6x1 vector)
F_e = zeros(6,1);

%section corresponding to x forces:
F_e(1) = h/3*f(1) + h/6*f(2);
F_e(2) = h/6*f(1) + h/3*f(2);

%section corresponding to y forces:
F_e(3) = 7*h/20*q(1) + 3*h/20*q(2);
F_e(4) = -h^2/20*q(1) - h^2/30*q(2);
F_e(5) = 3*h/20*q(1) + 7*h/20*q(2);
F_e(6) = h^2/30*q(1) + h^2/20*q(2);

%Rearrange vector to match structure of desired solution vector:
F_e = vect_rearrange( F_e );
