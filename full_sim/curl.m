function vflx = curl( stfn, stbc, parms, mats )

%Compute velocity flux from streamfunction.
%   Note: requires streamfunction from coarser grid on edge 
%         of current domain (stored in stbc) 

m = parms.m; n = parms.n;


%First get contribution without accounting for BCs:
vflx = mats.C * stfn;

%**Now for BCs

    %!!x-velocity block (terms on bottom and top edges):

        %Bottom part
            
            vflx( 1:m-1 ) = vflx( 1:m-1 ) - stbc( 1:m-1 ) ;

        %Top part

            %indices for vel grid
            topf = (n-1)*(m-1) + (1 : m-1);
            
            %indices for stfcn grid
            tops = (n-2)*(m-1) + (1 : m-1);
            
            %incorporate into curl
            vflx( topf ) = vflx( topf ) + stbc( tops );

    %!!    

    %!!y-velocity block

        %yvel index begins after x vels
        nadd = (m-1) * n; 
        
        %left part

            %indices for vel grid
            leftf = nadd + (1 : m : m*(n-2) + 1) ;
            
            %indices for strmfcn
            lefts = ( 1 : m-1 : (n-2)*(m-1) + 1 );

            %incorporate into curl
            vflx( leftf ) = vflx( leftf ) + stbc( lefts ) ;


        %right part

            %indices for vel grid
            rightf = nadd + ( m : m : m*(n-2) + m );
            
            %indices for stfcn grid
            rights = ( m-1 : m-1 : (n-1)*(m-1) );

            %points that need to average coarser domain:
            vflx( rightf ) = vflx( rightf ) - stbc( rights );

    %!!

%**

