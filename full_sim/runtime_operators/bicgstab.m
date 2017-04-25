function x = bicgstab( x, b, parms, mats )

%Solve Ax = b using the biCG-stabilized method

% display('Entering biCGstab...')

dt = parms.dt;

err = 1;
eps = 1e-14;

iter = 0;


B2 = 2/dt * mats.Itilde_flag * mats.sol_mat * mats.Q_flag * mats.W_flag;
h = parms.len / parms.m;
B2 = B2 * h / parms.ds; %scaling associated with ratio of IB to flow spacing
B2 = B2 / dt ; %Scaling associated with time stepping


r = b - b_times( x, parms, mats ) - B2 * x;
rhat = r;

rho_o = 1;
alpha = 1;
om = 1;

nu = 0;
p = 0;

eps = eps / h;

while ( err >= eps & iter < 1000 )
    
   rho_n = rhat' * r;
   bta = (rho_n / rho_o ) * (alpha / om );
   
   %Update rho
   rho_o = rho_n;
   
   p = r + bta * (p - om * nu );
   
   nu = b_times( p, parms, mats ) + B2 * p;
   
   alpha = rho_n / ( rhat' * nu );
   
   h_cg = x + alpha * p;
   
   sv = r - alpha * nu;
   
   tv = b_times( sv, parms, mats ) + B2 * sv;
   
   om = (tv' * sv) / (tv' * tv);
   
   x = h_cg + om * sv;
   
   r = sv - om * tv;
   
   err = r'*r;
   
   iter = iter + 1;
    
    
end


if (iter == 1000)
   
    warning('biCG terminated before convergence!')
end
