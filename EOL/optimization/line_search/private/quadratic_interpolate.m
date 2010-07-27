function [t, val] = quadratic_interpolate(t0, t1, f0, f1, d0)
h = t1 - t0;

% quadratic polinomial coefficients
p = [f1 - f0 - h*d0, d0, f0];

% location of minimimum
if f1 - f0 - h*d0 > 0 
    y = -0.5*d0*h/(f1 - f0 - h*d0);
else
    error('ooops');
end

% convert from y to t
t = y*h + t0;
val =  polyval(p, t);
