function W = get_W_main( W, glev, grid_parms )

%Build part of W corresponding to gridlevel glev (without accounting 
%for BC's). Note that grid_parms contains info on m, n, mg, and number of
%rows and columns of W.

m = grid_parms.m; n = grid_parms.n; mg = grid_parms.mg; 
len = grid_parms.len;
nrows = grid_parms.nrows;
ncols = grid_parms.ncols;

%y-velocity index starts at the end of all x-velocities
n_add = get_velx_ind( m-1, n, mg, grid_parms );

%grid spacing on finest grid level
del = len / m;

%grid spacing on current grid
delb = del * 2.d0^(glev-1);
    
    
%1st gridlevel is different from others:
if glev == 1
    
    %-- x-vel block
    
        %vorticity points above x-velocity point
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

        indvortx = repmat(1 : m-1,[1,n-1]);
        indvorty = repelem(1 : n-1, m-1);
        vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

        W = W + 1/2*sparse( velx_ind, vort_ind, ...
            ones(size(vort_ind)), nrows, ncols)/ (delb^2);
        
        %vorticity points below x-velocity point
        indvelx = repmat(1 : m-1, [1, n-1]);
        indvely = repelem(2 : n, m-1);
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

        indvortx = repmat(1 : m-1, [1, n-1]);
        indvorty = repelem(1 : n-1, m-1);
        vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

        W = W + 1/2*sparse( velx_ind, vort_ind,...
            ones(size(vort_ind)), nrows, ncols)/ (delb^2);

    %--
    
    %-- y-vel block
    
        %vorticity points to the right of y-velocity point
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

        indvortx = repmat(1 : m-1,[1,n-1]);
        indvorty = repelem(1 : n-1, m-1);
        vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

        W = W + 1/2*sparse(n_add + vely_ind, vort_ind, ...
            ones(size(vort_ind)), nrows, ncols)/ (delb^2);

        %vorticity points to the left of y-velocity point
        indvelx = repmat(2 : m,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

        indvortx = repmat(1 : m-1,[1,n-1]);
        indvorty = repelem(1 : n-1, m-1);
        vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

        W = W + 1/2*sparse(n_add + vely_ind, vort_ind, ...
            ones(size(vort_ind)), nrows, ncols)/ (delb^2);
    
    %--
 
%larger grid levels
else
    
    %-- x-vel blocks
        
        %vorticity points above x-velocity point
        
            %bottom region that doesn't overlap with finer grid
            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1,n/4]);
            indvorty = repelem(1 : n/4, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            W = W + 1/2*sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols)/ (delb^2);
            
            %top region that doesn't overlap with finer grid
            indvelx = repmat(1 : m-1,[1, (n/4-1)] );
            indvely = repelem( 3*n/4+1 : n-1, m-1 );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat( 1 : m-1,[1, (n/4-1)] );
            indvorty = repelem( 3*n/4+1 : n-1, m-1 );
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            W = W + 1/2*sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols) / (delb^2);
            
            %left and right parts that do have overlap
            indvelx = repmat([ 1 : m/4, 3*m/4 : m-1],[1, (n/2)] );
            indvely = repelem( n/4+1 : 3*n/4, m/2 );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat([ 1 : m/4, 3*m/4 : m-1],[1, (n/2)] );
            indvorty = repelem( n/4+1 : 3*n/4, m/2 );
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            W = W + 1/2*sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols)/ (delb^2);
            
            
        %vorticity points below x-velocity point
        
            %bottom region that doesn't overlap with finer grid
            indvelx = repmat(1 : m-1,[1,n/4-1]);
            indvely = repelem(2 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1,n/4-1]);
            indvorty = repelem(1 : n/4-1, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            W = W + 1/2*sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols)/ (delb^2);
            
            %top region that doesn't overlap with finer grid
            indvelx = repmat(1 : m-1,[1, n/4] );
            indvely = repelem( 3*n/4+1 : n, m-1 );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat( 1 : m-1,[1, n/4] );
            indvorty = repelem( 3*n/4 : n-1, m-1 );
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            W = W + 1/2*sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols)/ (delb^2);
            
            %left and right parts that do have overlap
            indvelx = repmat([ 1 : m/4, 3*m/4 : m-1],[1, n/2] );
            indvely = repelem( n/4+1 : 3*n/4, m/2 );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat([ 1 : m/4, 3*m/4 : m-1],[1, n/2] );
            indvorty = repelem( n/4 : 3*n/4-1, m/2 );
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            W = W + 1/2*sparse( velx_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols)/ (delb^2);
    
    
    %--
    
    %-- y-vel blocks
    
        %vorticity points to the right of y-velocity point
        
            %bottom part (no overlap with fine grid)
            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1,n/4]);
            indvorty = repelem(1 : n/4, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            W = W + 1/2*sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols)/ (delb^2);
            
            
            %top part (no overlap with fine grid)
            indvelx = repmat(1 : m-1,[1, n/4]);
            indvely = repelem(3*n/4 : n-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1, n/4]);
            indvorty = repelem(3*n/4 : n-1, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            W = W + 1/2*sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols)/ (delb^2);
            
            %middle part (overlap with fine grid)
            indvelx = repmat([1 : m/4, 3*m/4+1 : m-1],[1,3*n/4-(n/4+1)]);
            indvely = repelem(n/4+1 : 3*n/4-1, (m/4 + m/4-1) );
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat([1 : m/4, 3*m/4+1 : m-1],[1,3*n/4-(n/4+1)]);
            indvorty = repelem(n/4+1 : 3*n/4-1, (m/4 + m/4-1) );
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            W = W + 1/2*sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols)/ (delb^2);
            
            

        %vorticity points to the left of y-velocity point
        
            %bottom part (no overlap with fine grid)
            indvelx = repmat(2 : m,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1,n/4]);
            indvorty = repelem(1 : n/4, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            W = W + 1/2*sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols)/ (delb^2);
            
            %top part (no overlap with fine grid)
            indvelx = repmat(2 : m,[1,n/4]);
            indvely = repelem(3*n/4 : n-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat(1 : m-1,[1,n/4]);
            indvorty = repelem(3*n/4 : n-1, m-1);
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            W = W + 1/2*sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols)/ (delb^2);
            
            %middle part (overlap with fine grid)
            indvelx = repmat([ 2:m/4, 3*m/4+1:m],[1,n/2-1]);
            indvely = repelem( n/4+1 : 3*n/4-1, (m/4-1 + m/4) );
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvortx = repmat([1 : m/4-1, 3*m/4:m-1],[1,n/2-1]);
            indvorty = repelem( n/4+1 : 3*n/4-1, (m/4-1 + m/4));
            vort_ind = get_vort_ind(indvortx,indvorty,glev,grid_parms);

            W = W + 1/2*sparse(n_add + vely_ind, vort_ind, ...
                ones(size(vort_ind)), nrows, ncols)/ (delb^2);
    
    
    %--
    
    
    
end


