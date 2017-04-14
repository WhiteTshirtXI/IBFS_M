function [xb, ds] = build_cylinder( r, h )

%build cylinder of radius r using flow grid spacing of h
%(cylinder will have a spacing of 2h)

circum = 2 * pi * r; %Circumference of the circle


%Get # of points such that ds = 2h
n = floor( circum / h / 2 );

int =  2*pi / n ;
spt = 0 : int : (n-1)*int;
xfun = @(z) r*cos(z);
yfun = @(z) r*sin(z);
xhat = xfun(spt);
yhat = yfun(spt);

xb = [xhat, yhat];

save('body.mat', 'xb' )

%sanity check: make sure ds is equal to 2 * h
ds = sqrt( (xhat(2) - xhat(1))^2 + (yhat(2) - yhat(1))^2 ) ;
