function params = get_inds( params )

m = params.m; n = params.n;

%--Get indices for outer boundary of grid:
    %_x corresponds to x vel index and _y to y vel index
    
    bott_y = 1 : m;
    left_x = (m+1)*(0 : n-1) + 1;
    right_x = (m+1) + (m+1)*(0 : n-1);
    top_y = (m*n) + (1: m);
    
    bndryx = [left_x right_x];
    bndryy = [bott_y top_y];
    
    
%--

%--Indices that overlap with finer grid:
    %x vels first
    in_leftx = m/4 + 3;
    in_rightx = 3*m/4 - 1 ;
    in_downx = n/4 + 2;
    in_upx = 3*n/4 - 1;
    
    %y vels
    in_lefty = m/4 +2;
    in_righty = 3*m/4 - 1;
    in_downy = n/4 + 3;
    in_upy = 3*n/4 - 1;
    
    %Store them in inner...
    inx = []; iny = [];
    for j = in_downx : in_upx
        inx = [inx, (j-1)*(m+1) + (in_leftx : in_rightx)];
        
    end
    
    for j = in_downy : in_upy
        iny = [iny, (j-1)*m + (in_lefty : in_righty) ];
        
    end
        
params.bndryx = bndryx; params.bndryy = bndryy; params.inx = inx;
params.iny = iny;
    
%--

%Compute lengths of truncated velocity vectors:

params.nmg1 = length(params.indbgx)-length(params.bndryx) + ...
    (length(params.indbgy)- length(params.bndryy));

params.nmid = length(params.indbgx) - (length(params.bndryx) + length(params.inx)) + ...
        length(params.indbgy) - (length(params.bndryy) + length(params.iny));

params.nmglst = length(params.indbgx) - length(params.inx) + ...
    length(params.indbgy) - length(params.iny);


%Total number of points in velocity vectors:
params.nvel_tot = params.nmg1 + (params.mg - 1) * params.nmid; 


