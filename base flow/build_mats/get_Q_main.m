function Q = get_Q_main( Q, glev, grid_parms )

%Build part of Q corresponding to gridlevel glev (without accounting 
%for BC's). Note that grid_parms contains info on m, n, mg, and number of
%rows and columns of Q.

m = grid_parms.m; n = grid_parms.n; mg = grid_parms.mg;
nrows = grid_parms.nrows;
ncols = grid_parms.ncols;

%y-velocity index starts at the end of all x-velocities
n_add = get_velx_ind( m-1, n, mg, grid_parms );


%1st gridlevel is different from others:
if glev == 1
    
    %--Block corresponding to x-velocities
    
        %y-vels to bottom left of current x-vel point
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(2 : n, m-1);
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);
        
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);
        
        Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
            ones(size(vely_ind)), nrows, ncols);
        
        %y-vels to bottom right of current x-vel point
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(2 : n, m-1);
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);
        
        indvelx = repmat(2 : m,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);
        
        Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
            ones(size(vely_ind)), nrows, ncols);
        
        %y-vels to top left of current x-vel point
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);
        
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);
        
        Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
            ones(size(vely_ind)), nrows, ncols);
        
        %y-vels to top right of current x-vel point
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);
        
        indvelx = repmat(2 : m,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);
        
        Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
            ones(size(vely_ind)), nrows, ncols);
    
    %--
    
    %--Block corresponding to y-velocities
    
        %x-vels to bottom left of current y-vel point
        indvelx = repmat(2 : m,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);
        
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);
        
        Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
            ones(size(vely_ind)), nrows, ncols);
        
        %x-vels to bottom right of current y-vel point
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);
        
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);
        
        Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
            ones(size(vely_ind)), nrows, ncols);
        
        %x-vels to top left of current y-vel point
        indvelx = repmat(2 : m,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);
        
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(2 : n, m-1);
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);
        
        Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
            ones(size(vely_ind)), nrows, ncols);
        
        %x-vels to top right of current y-vel point
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(1 : n-1, m-1);
        vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);
        
        indvelx = repmat(1 : m-1,[1,n-1]);
        indvely = repelem(2 : n, m-1);
        velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);
        
        Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
            ones(size(vely_ind)), nrows, ncols);
    
    %--
    
else
    
    %--Block corresponding to x-velocities
    
        %y-vels to bottom left of current x-vel point
        
            %bottom part (no overlap with fine grid)
            indvelx = repmat(1 : m-1,[1,n/4-1]);
            indvely = repelem(2 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(1 : m-1,[1,n/4-1]);
            indvely = repelem(1 : n/4-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %top part (no overlap with fine grid)
            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(3*n/4+1 : n, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(3*n/4 : n-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %middle part (overlap with fine grid)
            indvelx = repmat([1:m/4, 3*m/4+1:m-1],[1,n/2]);
            indvely = repelem(n/4+1 : 3*n/4, (m/4+m/4-1));
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat([1:m/4, 3*m/4+1:m-1],[1,n/2]);
            indvely = repelem(n/4 : 3*n/4-1, (m/4+m/4-1));
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %add bottom right corner of fine grid
            indvelx = 3*m/4;
            indvely = n/4+1;
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = 3*m/4;
            indvely = n/4;
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
        %y-vels to bottom right of current x-vel point
        
            %bottom part (no overlap with fine grid)
            indvelx = repmat(1 : m-1,[1,n/4-1]);
            indvely = repelem(2 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(2 : m,[1,n/4-1]);
            indvely = repelem(1 : n/4-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %top part (no overlap with fine grid)
            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(3*n/4+1 : n, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(2 : m,[1,n/4]);
            indvely = repelem(3*n/4 : n-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %middle part (overlap with fine grid)
            indvelx = repmat([1:m/4-1, 3*m/4:m-1],[1,n/2]);
            indvely = repelem(n/4+1 : 3*n/4, (m/4-1 + m/4));
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat([2:m/4, 3*m/4+1:m],[1,n/2]);
            indvely = repelem(n/4 : 3*n/4-1, (m/4-1 + m/4));
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %bottom left corner of fine grid
            indvelx = m/4;
            indvely = n/4+1;
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = m/4+1;
            indvely = n/4;
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
        %y-vels to top right of current x-vel point
        
            %bottom part (no overlap with fine grid)
            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(2 : m,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %top part (no overlap with fine grid)
            indvelx = repmat(1 : m-1,[1,n/4-1]);
            indvely = repelem(3*n/4+1 : n-1, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(2 : m,[1,n/4-1]);
            indvely = repelem(3*n/4+1 : n-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %middle part (overlap with fine grid)
            indvelx = repmat([1:m/4-1, 3*m/4:m-1],[1,n/2]);
            indvely = repelem(n/4+1 : 3*n/4, (m/4-1 + m/4));
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat([2:m/4, 3*m/4+1:m],[1,n/2]);
            indvely = repelem(n/4 +1: 3*n/4, (m/4-1 + m/4));
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %top left corner of fine grid
            indvelx = m/4;
            indvely = 3*n/4;
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = m/4+1;
            indvely = 3*n/4;
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
        %y-vels to top left of current x-vel point
        
            %bottom part (no overlap with fine grid)
            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %top part (no overlap with fine grid)
            indvelx = repmat(1 : m-1,[1,n/4-1]);
            indvely = repelem(3*n/4+1 : n-1, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(1 : m-1,[1,n/4-1]);
            indvely = repelem(3*n/4+1 : n-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %middle part (overlap with fine grid)
            indvelx = repmat([1:m/4, 3*m/4+1:m-1],[1,n/2]);
            indvely = repelem(n/4+1 : 3*n/4, (m/4 + m/4-1 ));
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat([1:m/4, 3*m/4+1:m-1],[1,n/2]);
            indvely = repelem(n/4 +1: 3*n/4, (m/4 + m/4-1 ));
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %top left corner of fine grid
            indvelx = 3*m/4;
            indvely = 3*n/4;
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            indvelx = 3*m/4;
            indvely = 3*n/4;
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            Q = Q - 1/4*sparse( velx_ind, n_add + vely_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
     %--
     
     %--Blocks corresponding to y-velocities
     
        %x-vels to bottom left of current y-vel point
        
            %bottom (no overlap with finer grid)
            indvelx = repmat(2 : m,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %top (no overlap with finer grid)
            indvelx = repmat(2 : m,[1,n/4-1]);
            indvely = repelem(3*n/4 +1 : n-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(1 : m-1,[1,n/4-1]);
            indvely = repelem(3*n/4 +1 : n-1, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %middle (overlap with finer grid)
            indvelx = repmat([2 : m/4, 3*m/4+1:m],[1,n/2]);
            indvely = repelem(n/4 +1 : 3*n/4, (m/4-1 +m/4) );
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat([1 : m/4-1, 3*m/4:m-1],[1,n/2]);
            indvely = repelem(n/4 +1 : 3*n/4, (m/4-1 +m/4) );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);
        
            %top left corner of finer grid
            indvelx = m/4+1;
            indvely = 3*n/4;
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = m/4;
            indvely = 3*n/4;
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);
 
        %x-vels to bottom right of current y-vel point
        
            %bottom (no overlap with finer grid)
            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %top (no overlap with finer grid)
            indvelx = repmat(1 : m-1,[1,n/4-1]);
            indvely = repelem(3*n/4 +1 : n-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(1 : m-1,[1,n/4-1]);
            indvely = repelem(3*n/4 +1 : n-1, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %middle (overlap with finer grid)
            indvelx = repmat([1 : m/4, 3*m/4+1:m-1],[1,n/2]);
            indvely = repelem(n/4 +1 : 3*n/4, (m/4 +m/4-1) );
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat([1 : m/4, 3*m/4+1:m-1],[1,n/2]);
            indvely = repelem(n/4 +1 : 3*n/4, (m/4 + m/4-1 ) );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);
        
            %top right corner of finer grid
            indvelx = 3*m/4;
            indvely = 3*n/4;
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = 3*m/4;
            indvely = 3*n/4;
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);           
            
            
        %x-vels to top right of current y-vel point
        
            %bottom (no overlap with finer grid)
            indvelx = repmat(1 : m-1,[1,n/4-1]);
            indvely = repelem(1 : n/4-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(1 : m-1,[1,n/4-1]);
            indvely = repelem(2 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %top (no overlap with finer grid)
            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(3*n/4 : n-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(3*n/4 +1 : n, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %middle (overlap with finer grid)
            indvelx = repmat([1 : m/4, 3*m/4+1:m-1],[1,n/2]);
            indvely = repelem(n/4 : 3*n/4-1, (m/4 +m/4-1) );
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat([1 : m/4, 3*m/4+1:m-1],[1,n/2]);
            indvely = repelem(n/4 +1 : 3*n/4, (m/4 + m/4-1 ) );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);
        
            %bottom right corner of finer grid
            indvelx = 3*m/4;
            indvely = n/4;
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = 3*m/4;
            indvely = n/4+1;
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);  
            
        %x-vels to top left of current y-vel point
        
            %bottom (no overlap with finer grid)
            indvelx = repmat(2 : m,[1,n/4-1]);
            indvely = repelem(1 : n/4-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(1 : m-1,[1,n/4-1]);
            indvely = repelem(2 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %top (no overlap with finer grid)
            indvelx = repmat(2 : m,[1,n/4]);
            indvely = repelem(3*n/4 : n-1, m-1);
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(3*n/4 +1 : n, m-1);
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);
            
            %middle (overlap with finer grid)
            indvelx = repmat([2 : m/4, 3*m/4+1:m],[1,n/2]);
            indvely = repelem(n/4 : 3*n/4-1, (m/4 +m/4-1) );
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = repmat([1 : m/4-1, 3*m/4:m-1],[1,n/2]);
            indvely = repelem(n/4 +1 : 3*n/4, (m/4 + m/4-1 ) );
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);
        
            %bottom left corner of finer grid
            indvelx = m/4+1;
            indvely = n/4;
            vely_ind = get_vely_ind(indvelx,indvely,glev,grid_parms);

            indvelx = m/4;
            indvely = n/4+1;
            velx_ind = get_velx_ind(indvelx,indvely,glev,grid_parms);

            Q = Q + 1/4*sparse( n_add + vely_ind, velx_ind, ...
                ones(size(vely_ind)), nrows, ncols);           

     %--
    
end


