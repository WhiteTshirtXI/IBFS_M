function [q, s] = circ2_st_vflx( gamma, ngrids, parms, mats )

%Take vorticity on ngrids and return vel flux and streamfunction on those
%grids


m = parms.m; n = parms.n; 
nq = (m-1)*n + (n-1) * m; ngam = (m-1)*(n-1);

q = zeros( nq, ngrids );
s = zeros( ngam, ngrids );

for j = ngrids : -1 : 1

    %--Solve Poisson problem for stfn
    
        %BCs for Poisson problem (will be overwritten if mg > 1)
        stbc = zeros( ngam, 1 );
    
        %don't need bcs for largest grid
        if ( j == ngrids ) 
            
            %Solve for streamfcn
            s(:,j) = RCinv( gamma(:,j), parms, mats );
            
        %get bcs for other grid levels
        else
            
            %Get bcs from streamfunction on coarser gridlevel
            stbc = get_stfn_BCs( stbc, s(:,j+1), parms );
            
            %Solve for streamfcn
            s(:,j) = RCinv( gamma(:,j) + stbc, parms, mats );
            
        end

    %--

    %--Get velocity on first grid from stream function

        q = curl( stfn, stbc, parms );

    %--
    
end

