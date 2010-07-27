function [t, val] = cubic_interpolate(t0,t1, f0, f1, d0, d1)
h = t1 - t0;

% cubic polinomial coefficients
% a*y^3 + b*y^2 + c*y + d
p = genCubicPolynomialCoeff([f0, f1], [d0*h, d1*h]);
y = findCubicPolyMinimum(p);

t = y*h + t0;
val = polyval(p, t);
