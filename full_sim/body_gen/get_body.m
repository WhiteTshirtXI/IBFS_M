function [parms, soln] = get_body( parms )

%Build body for simulation from one of a number of labels

%'cyl' = cylinder with radius 1/2
%'flg' = flag with length 1

 if parms.body == 'cyl'
        
    r = 0.5;
    
    [soln.xb, parms.ds] = build_cylinder( r, parms.len/parms.m );
    
    parms.nb = length(soln.xb) / 2;
    
     
        
 elseif parms.body == 'flg'
     
    h = parms.len / parms.m;
    parms.ds = 2*h;

    %Reference configuration for flag
    xb0x = 0 : parms.ds : 1;
    xb0y = zeros( size( xb0x ) );

    parms.xb0 = [xb0x xb0y];
    parms.nb = length( xb0x );
    soln.xb = parms.xb0; %body position begins in undeformed state
    
 else
     
     if ( exist('body.mat') ~= 2 )
     
         error(['Body type not supported! Have to supply your own body'],...
            ['and store it in body.mat with the variables xb, ds, xb0, and nb'])
        
     else
         
         load('body.mat')
         
         soln.xb = xb;
         parms.ds = ds;
         parms.nb = nb;
         
         if (exist( 'xb0') == 1 )
             parms.xb0 = xb0;
         end
    
     end
     
 end
 



