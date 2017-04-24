function index = get_ind_flag( iel )

nnel = 2; % # of nodes per element
ndof = 3; % # of DOF per node

edof = nnel * ndof; %degrees of freedom for an element

start = (iel - 1) * (nnel - 1) * ndof; 

index = zeros( edof, 1);
for i = 1 : edof
   
    index(i) = start + i;
    
end