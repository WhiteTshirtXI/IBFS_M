function qtr = truncate( q, params )


%Take velocity matrix of size nq x mg as given by fortran code and truncate
%to eliminate redundant information on overlapping grids. Also, convert the
%velocity as given by the code into physical units.

nq = params.nq; mg = params.mg ;  len = params.len;
m = params.m; n = params.n; inx = params.inx; iny = params.iny;
bndryx = params.bndryx; bndryy = params.bndryy; 
indbgx = params.indbgx; indbgy = params.indbgy;

%Grid 1:
qbg = q(1 : nq );

% %Scale to get velocities from fluxes:
% fac = 1;
% dlt = len ./ m *fac;
% qbg = qbg ./ dlt;

qx = qbg( 1 : (m+1)*n); qy = qbg( (m+1)*n + 1 : nq );


indx = indbgx;
indy = indbgy;
indx(bndryx) = 0;
indy(bndryy) = 0;
qtrx = qx( indx ~= 0);
qtry = qy( indy ~= 0 );
qtr = [qtrx; qtry];

if mg > 1
    
    for j = 2 : mg
        
        qbg = q( (j-1)*nq + (1 : nq) );
        
%         %Scale to get velocities from fluxes:
%         fac = 2^(j-1);
%         dlt = len ./ m *fac;
%         qbg = qbg ./ dlt;

        qx = qbg( 1 : (m+1)*n); qy = qbg( (m+1)*n + 1 : nq );
        
        indx = indbgx;
        indx(bndryx) = 0;
        indx(inx) = 0;
        qtrxtmp = qx(indx ~= 0 );
        
        indy = indbgy;
        indy(bndryy) = 0;
        indy(iny) = 0;
        qtrytmp = qy(indy ~= 0 );
        
        qtr = [qtr; qtrxtmp; qtrytmp];
        
        
    end
    
        
%     %Last gridlevel:
%     qbg = q( (mg-1)*nq + (1 : nq) );
%     
%     %Scale to get velocities from fluxes:
%     fac = 2^(mg-1);
%     dlt = len ./ m *fac;
%     qbg = qbg ./ dlt;
%         
%     qx = qbg( 1 : (m+1)*n); qy = qbg( (m+1)*n + 1 : nq );
% 
%     indx = indbgx;
%     indx(inx) = 0;
%     qtrxtmp = qx( indx ~= 0 );
%     
%     indy = indbgx;
%     indy(iny) = 0;
%     qtrytmp = qy( indy ~= 0 );
%     qtr = [qtr; qtrxtmp; qtrytmp];
    
    
    
end





