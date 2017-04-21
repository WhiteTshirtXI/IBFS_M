function x = b_times( z, parms, mats )

%Performs one matrix multiply of B*z, where B is the matrix used to solve
%for the surface stresses that enforce the no-slip boundary condition.

% (B arises from an LU factorization of the full system)

%--Initialize

    m = parms.m; n = parms.n; mg = parms.mg;
    h = (parms.len / parms.m); ngam = (m-1)* (n-1);

    circ = zeros( ngam, mg ); 
    
%--

%-- get circulation from surface stress

    %Get circ on 1st grid level
    circ(:, 1) = Ainv( mats.R * mats.ET * z, 1, parms, mats );
    %We don't include BCs from coarse grid for Ainv because ET*z is compact

    %Coarsify circulation to second grid level to get BCs for stfn
    if mg > 1
        circ(:,2) = coarsify( circ(:,1), circ(:,2), parms );
    end
%--

%-- get vel flux from circulation

    %only need to work with 2 grid levels
    if parms.mg > 1
        vflx = circ2_st_vflx( circ, 2,  parms, mats);
    else
        vflx = circ2_st_vlx( circ, 1, parms, mats );
    end

%--
    
%--Interpolate onto the body and scale by h
    
    x = mats.E * vflx(:,1) / h;
    
%--    
