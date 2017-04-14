function C = get_C_main( C, glev, grid_parms )

m = grid_parms.m; n = grid_parms.n; mg = grid_parms.mg;
nrows = grid_parms.nrows;
ncols = grid_parms.ncols;

%1st gridlevel is different from the others
if glev == 1
    
    %--First build block corresponding to x-velocities

        %vorticity points above x-velocity point
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

        indvortx = repmat(1 : m-1,[1,n-1]);
        indvorty = repelem(1 : n-1, m-1);
        vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

        C = C + sparse( velx_ind, vort_ind, ...
            ones(size(vort_ind)), nrows, ncols);

        %vorticity points below x-velocity point
        indvelx = repmat(1 : m-1, [1, n-1]);
        indvely = repelem(2 : n, m-1);
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

        indvortx = repmat(1 : m-1, [1, n-1]);
        indvorty = repelem(1 : n-1, m-1);
        vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

        C = C - sparse( velx_ind, vort_ind,...
            ones(size(vort_ind)), nrows, ncols);

    %--    

    %--Now build y-velocity block

        %rows start at end of x-velocity block:
        n_add = get_velx_ind( m-1, n, mg, grid_parms );

        %vorticity points to the right of y-velocity point
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

        indvortx = repmat(1 : m-1,[1,n-1]);
        indvorty = repelem(1 : n-1, m-1);
        vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

        C = C - sparse(n_add + vely_ind, vort_ind, ...
            ones(size(vort_ind)), nrows, ncols);

        %vorticity points to the left of y-velocity point
        indvelx = repmat(2 : m,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

        indvortx = repmat(1 : m-1,[1,n-1]);
        indvorty = repelem(1 : n-1, m-1);
        vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

        C = C + sparse(n_add + vely_ind, vort_ind, ...
            ones(size(vort_ind)), nrows, ncols);

    %--
 
%larger grid levels    
elseif glev < mg
    
    %--First build block corresponding to x-velocities
    
        %vorticity points above x-velocity point
        
            %bottom region that doesn't overlap with finer grid
            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1,n/4]);
            indvorty = repelem(1 : n/4, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C + sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            %top region that doesn't overlap with finer grid
            indvelx = repmat(1 : m-1,[1, (n-1)-(3*n/4)] );
            indvely = repelem( 3*n/4+1 : n-1, m-1 );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat( 1 : m-1,[1, (n-1)-(3*n/4)] );
            indvorty = repelem( 3*n/4+1 : n-1, m-1 );
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C + sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            %left and right parts that do have overlap
            indvelx = repmat([ 1 : m/4, 3*m/4 : m-1],[1, (3*n/4-n/4)] );
            indvely = repelem( n/4+1 : 3*n/4, m/2 );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat([ 1 : m/4, 3*m/4 : m-1],[1, (3*n/4-n/4)] );
            indvorty = repelem( n/4+1 : 3*n/4, m/2 );
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C + sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
    
            
        %vorticity points below x-velocity point
        
            %bottom region that doesn't overlap with finer grid
            indvelx = repmat(1 : m-1,[1,n/4-1]);
            indvely = repelem(2 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1,n/4-1]);
            indvorty = repelem(1 : n/4-1, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C - sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            %top region that doesn't overlap with finer grid
            indvelx = repmat(1 : m-1,[1, n-(3*n/4)] );
            indvely = repelem( 3*n/4+1 : n, m-1 );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat( 1 : m-1,[1, n-(3*n/4)] );
            indvorty = repelem( 3*n/4 : n-1, m-1 );
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C - sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            %left and right parts that do have overlap
            indvelx = repmat([ 1 : m/4, 3*m/4 : m-1],[1, (3*n/4-n/4)] );
            indvely = repelem( n/4+1 : 3*n/4, m/2 );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat([ 1 : m/4, 3*m/4 : m-1],[1, (3*n/4-n/4)] );
            indvorty = repelem( n/4 : 3*n/4-1, m/2 );
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C - sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            
            
    %--Now build y-velocity block

        %rows start at end of x-velocity block:
        n_add = get_velx_ind( m-1, n, mg, grid_parms );

        %vorticity points to the right of y-velocity point
        
            %bottom part (no overlap with fine grid)
            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1,n/4]);
            indvorty = repelem(1 : n/4, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C - sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            
            %top part (no overlap with fine grid)
            indvelx = repmat(1 : m-1,[1, n/4]);
            indvely = repelem(3*n/4 : n-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1, n/4]);
            indvorty = repelem(3*n/4 : n-1, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C - sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            %middle part (overlap with fine grid)
            indvelx = repmat([1 : m/4, 3*m/4+1 : m-1],[1,3*n/4-(n/4+1)]);
            indvely = repelem(n/4+1 : 3*n/4-1, (m/4 + m/4-1) );
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat([1 : m/4, 3*m/4+1 : m-1],[1,3*n/4-(n/4+1)]);
            indvorty = repelem(n/4+1 : 3*n/4-1, (m/4 + m/4-1) );
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C - sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            

        %vorticity points to the left of y-velocity point
        
            %bottom part (no overlap with fine grid)
            indvelx = repmat(2 : m,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1,n/4]);
            indvorty = repelem(1 : n/4, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C + sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            %top part (no overlap with fine grid)
            indvelx = repmat(2 : m,[1,n/4]);
            indvely = repelem(3*n/4 : n-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1,n/4]);
            indvorty = repelem(3*n/4 : n-1, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C + sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            %middle part (overlap with fine grid)
            indvelx = repmat([ 2:m/4, 3*m/4+1:m],[1,n/2-1]);
            indvely = repelem( n/4+1 : 3*n/4-1, (m/4-1 + m/4) );
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat([1 : m/4-1, 3*m/4:m-1],[1,n/2-1]);
            indvorty = repelem( n/4+1 : 3*n/4-1, (m/4-1 + m/4));
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C + sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
    %--
 
%largest domain is slightly different because of outflow boundary
else
    
       
    %--First build block corresponding to x-velocities
    
        %vorticity points above x-velocity point
        
            %bottom region that doesn't overlap with finer grid
            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1,n/4]);
            indvorty = repelem(1 : n/4, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C + sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            %top region that doesn't overlap with finer grid
            indvelx = repmat(1 : m-1,[1, (n-1)-(3*n/4)] );
            indvely = repelem( 3*n/4+1 : n-1, m-1 );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat( 1 : m-1,[1, (n-1)-(3*n/4)] );
            indvorty = repelem( 3*n/4+1 : n-1, m-1 );
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C + sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            %left and right parts that do have overlap
            indvelx = repmat([ 1 : m/4, 3*m/4 : m-1],[1, (3*n/4-n/4)] );
            indvely = repelem( n/4+1 : 3*n/4, m/2 );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat([ 1 : m/4, 3*m/4 : m-1],[1, (3*n/4-n/4)] );
            indvorty = repelem( n/4+1 : 3*n/4, m/2 );
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C + sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
    
            
        %vorticity points below x-velocity point
        
            %bottom region that doesn't overlap with finer grid
            indvelx = repmat(1 : m-1,[1,n/4-1]);
            indvely = repelem(2 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1,n/4-1]);
            indvorty = repelem(1 : n/4-1, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C - sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            %top region that doesn't overlap with finer grid
            indvelx = repmat(1 : m-1,[1, n-(3*n/4)] );
            indvely = repelem( 3*n/4+1 : n, m-1 );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat( 1 : m-1,[1, n-(3*n/4)] );
            indvorty = repelem( 3*n/4 : n-1, m-1 );
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C - sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            %left and right parts that do have overlap
            indvelx = repmat([ 1 : m/4, 3*m/4 : m-1],[1, (3*n/4-n/4)] );
            indvely = repelem( n/4+1 : 3*n/4, m/2 );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat([ 1 : m/4, 3*m/4 : m-1],[1, (3*n/4-n/4)] );
            indvorty = repelem( n/4 : 3*n/4-1, m/2 );
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C - sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            
            
    %--Now build y-velocity block

        %rows start at end of x-velocity block:
        n_add = get_velx_ind( m-1, n, mg, grid_parms );

        %vorticity points to the right of y-velocity point
        
            %bottom part (no overlap with fine grid)
            indvelx = repmat(1 : m,[1,n/4]);
            indvely = repelem(1 : n/4, m);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m,[1,n/4]);
            indvorty = repelem(1 : n/4, m);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C - sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            
            %top part (no overlap with fine grid)
            indvelx = repmat(1 : m,[1, n/4]);
            indvely = repelem(3*n/4 : n-1, m);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m,[1, n/4]);
            indvorty = repelem(3*n/4 : n-1, m);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C - sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            %middle part (overlap with fine grid)
            indvelx = repmat([1 : m/4, 3*m/4+1 : m],[1,3*n/4-(n/4+1)]);
            indvely = repelem(n/4+1 : 3*n/4-1, (m/4 + m/4) );
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat([1 : m/4, 3*m/4+1 : m],[1,3*n/4-(n/4+1)]);
            indvorty = repelem(n/4+1 : 3*n/4-1, (m/4 + m/4) );
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C - sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            

        %vorticity points to the left of y-velocity point
        
            %bottom part (no overlap with fine grid)
            indvelx = repmat(2 : m,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1,n/4]);
            indvorty = repelem(1 : n/4, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C + sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            %top part (no overlap with fine grid)
            indvelx = repmat(2 : m,[1,n/4]);
            indvely = repelem(3*n/4 : n-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1,n/4]);
            indvorty = repelem(3*n/4 : n-1, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C + sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
            %middle part (overlap with fine grid)
            indvelx = repmat([ 2:m/4, 3*m/4+1:m],[1,n/2-1]);
            indvely = repelem( n/4+1 : 3*n/4-1, (m/4-1 + m/4) );
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat([1 : m/4-1, 3*m/4:m-1],[1,n/2-1]);
            indvorty = repelem( n/4+1 : 3*n/4-1, (m/4-1 + m/4));
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            C = C + sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols);
            
    %--
    
end

