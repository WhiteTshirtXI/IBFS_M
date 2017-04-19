function x = b_times( z, parms, mats )

%Performs one matrix multiply of B*z, where B is the matrix used to solve
%for the surface stresses that enforce the no-slip boundary condition.

% (B arises from an LU factorization of the full system)


%--Initialize
    m = parms.m; n = parms.n; mg = parms.mg;
    nq = (m-1)*n + (n-1)*m; lev = 1; %body is in first grid level

    circ = zeros( (m-1)*(n-1), mg ); stfn = circ; vflx = zeros(nq, 1);
%--


%-- get circulation from surface stress

    %Get circ on 1st grid level
    circ(:, 1) = Ainv( mats.R * mats.ET * z, lev, parms, mats );
    %We don't include BCs from coarse grid for Ainv because ET*z is compact

    %Coarsify circulation to second grid level to get BCs for stfn
    if mg > 1
        circ(:,2) = coarsify( circ(:,1) );
    end
%--

%--Solve Poisson problem for stfn

    %Solve on coarser grid level to get bcs for first grid level
    if mg > 1

        stfn = RCinv( circ(:,2), parms, mats );

        %First need bcs from 
        stbc = zeros( ngam, 1 );

        %from explicit treatment of circulation
        stbc = get_stfn_BCs( stbc, stfn, 1, parms );
        % (1 because we are getting bcs for first grid level)     
    end
%--
    
    




