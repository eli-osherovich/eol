function t = quadratic_interpolate(t0, t1, f0, f1, d0)


% Copyright 2010 Eli Osherovich.


% Convert the original interval [t0, t1] to [0,1].
h = t1 - t0;
d0 = d0*h;

% Find quadratic polinomial coefficients: a*x^2 + b*x + c
a = f1 - f0 - d0;
b = d0;
%c = f0;

% Find minimum of the quadratic
if a > 0 
    y = -b/2/a;
else
    error('EOL:QUADRATIC_INTERPOLATE:NonConvexParabola', ...
        'The resulting parabola is not convex: no minimimum exists');
end

% Convert back to the original interval [t0, t1].
t = y*h + t0;
