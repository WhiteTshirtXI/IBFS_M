function circc = coarsify( circf, circc, parms )

%Take circulation on a finer grid (circf) and use it to replace the
%circulation on a coarser grid (circc) in the overlapping region.

m = parms.m; n = parms.n; 

%--Define indices of overlapping region for coarse grid:

    %1st row
    ind1 = n/4*(m-1) + (m/4 + 1 : 3*m/4 -1); 
    Ind1 = repmat( ind1, [1,n/2 - 1]);

    %indices to add for remaining rows
    indadd = 0 : m-1 : (n/2-2)*(m-1); 
    Indadd = repelem( indadd, length(ind1) );


    coarse_c = Ind1 + Indadd;

%--

%--Indices for various pieces of fine grid that will be used

    %!!center point
    
        %1st row
        ind1 = m-1 + ( 2 : 2 : m-2);
        Ind1 = repmat( ind1, [1, n/2 - 1] );
        
        %indices to add
        indadd = repelem( 0 : 2*(m-1) : 2*(n/2-2) *(m-1), length(ind1) );
        
        fine_c = Ind1 + indadd;
        
    
    %!!
    
    %!!top point
    
        %1st row
        ind1 = 2*(m-1) + ( 2 : 2 : m-2);
        Ind1 = repmat( ind1, [1, n/2 - 1] );
        
        %indices to add
        indadd = repelem( 0 : 2*(m-1) : 2*(n/2-2) *(m-1), length(ind1) );
        
        fine_t = Ind1 + indadd;
        
    
    %!!
    
    %!!bottom point
    
        %1st row
        ind1 = ( 2 : 2 : m-2);
        Ind1 = repmat( ind1, [1, n/2 - 1] );
        
        %indices to add
        indadd = repelem( 0 : 2*(m-1) : 2*(n/2-2) *(m-1), length(ind1) );
        
        fine_b = Ind1 + indadd;
        
    
    %!!

    
    %!!left point
        %1st row
        ind1 = (m-1) + ( 1 : 2 : m-3);
        Ind1 = repmat( ind1, [1, n/2 - 1] );

        %indices to add
        indadd = repelem( 0 : 2*(m-1) : 2*(n/2-2) *(m-1), length(ind1) );

        fine_l = Ind1 + indadd;
        
    %!!
    
    %!!right point
        %1st row
        ind1 = (m-1) + ( 3 : 2 : m-1);
        Ind1 = repmat( ind1, [1, n/2 - 1] );

        %indices to add
        indadd = repelem( 0 : 2*(m-1) : 2*(n/2-2) *(m-1), length(ind1) );

        fine_r = Ind1 + indadd;
        
    %!!
    
    %!!top left point
        %1st row
        ind1 = 2*(m-1) + ( 1 : 2 : m-3);
        Ind1 = repmat( ind1, [1, n/2 - 1] );

        %indices to add
        indadd = repelem( 0 : 2*(m-1) : 2*(n/2-2) *(m-1), length(ind1) );

        fine_tl = Ind1 + indadd;
        
    %!!
    
    %!!top right point
        %1st row
        ind1 = 2*(m-1) + ( 3 : 2 : m-1);
        Ind1 = repmat( ind1, [1, n/2 - 1] );

        %indices to add
        indadd = repelem( 0 : 2*(m-1) : 2*(n/2-2) *(m-1), length(ind1) );

        fine_tr = Ind1 + indadd;
        
    %!!

     %!!bottom left point
        %1st row
        ind1 = ( 1 : 2 : m-3);
        Ind1 = repmat( ind1, [1, n/2 - 1] );

        %indices to add
        indadd = repelem( 0 : 2*(m-1) : 2*(n/2-2) *(m-1), length(ind1) );

        fine_bl = Ind1 + indadd;
        
    %!!
    
    %!!bottom right point
        %1st row
        ind1 = ( 3 : 2 : m-1);
        Ind1 = repmat( ind1, [1, n/2 - 1] );

        %indices to add
        indadd = repelem( 0 : 2*(m-1) : 2*(n/2-2) *(m-1), length(ind1) );

        fine_br = Ind1 + indadd;
        
    %!!
    

%--

%--interpolate finer circulation onto coarser one

    circc( coarse_c ) = circf( fine_c ) + ...
        0.5 * ( circf( fine_t ) + circf( fine_b ) + ...
                circf( fine_l ) + circf( fine_r ) ) + ...
        0.25 * ( circf( fine_tl ) + circf( fine_tr ) + ...
                 circf( fine_bl ) + circf( fine_br ) ) ;

%--






