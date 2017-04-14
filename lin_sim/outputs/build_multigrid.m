clear all, close all, clc


%%
m = 12;
n = 8;
len = 2;

offsetx = 0.2;
offsety = 0.5;

clr = {'bx', 'r+', 'ko'};

for lev = 1 : 2

    fac = 2^(lev-1); 
    
    % Grid spacing in both directions for current grid
    delta = len ./ m *fac;              
    
    % Offset in x direction for current grid
    offx = 2^(lev-1) * len/2 - len/2 + offsetx ;    
    
    % Offset in y direction for current grid
    offy = 2^(lev-1) * (n*len/m)/2 - (n*len/m)/2 + offsety ;
    
    xv = 0 : delta : (m) * delta;
    yv = 0 : delta : (n) * delta;
    
    xv = xv - offx;
    yv = yv - offy;
    
    
    figure(10), hold on
    
    for jj = 1 : length( xv )
   
        plot( xv(jj) * ones(length(yv)), yv, clr{lev} )
        
    end
    
    
end

figure(10), set(gca,'Xtick',[],'Ytick',[],'Xticklabels',[],'Yticklabels',[])
axis([min(xv)-1 max(xv)+1 min(yv)-1 max(yv) + 1])
box on
