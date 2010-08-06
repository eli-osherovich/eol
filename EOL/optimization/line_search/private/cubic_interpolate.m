function t = cubic_interpolate(t0, t1, f0, f1, d0, d1)


% Copyright 2010 Eli Osherovich.



% Convert the original interval [t0, t1] to [0,1].
h = t1 - t0;
d0 = d0*h;
d1 = d1*h;

% Find cubic polinomial coefficients: a*y^3 + b*y^2 + c*y + d
% Using divided differences.
%d = f0;
c = d0;
b = 3*(f1 - f0) - 2*d0 - d1;
a = 2*(f0 - f1) + d0 + d1;

% Find minimimum of the cubic
% assuming a /= 0 !!!
if b < 0 % no canellation
    y = (sqrt(b^2-3*a*c) - b)/(3*a);
else % avoid cancellation
    y = -c/(sqrt(b^2 - 3*a*c) + b);
end
    
% Convert back to the original interval [t0, t1].
t = y*h + t0;
