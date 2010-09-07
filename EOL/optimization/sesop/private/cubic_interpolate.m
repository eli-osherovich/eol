function newt = cubic_interpolate(t0,t1, f0, f1, d0, d1)
h = t1 - t0;

% cubic polinomial coefficients
% a*y^3 + b*y^2 + c*y + d
p = genCubicPolynomialCoeff([f0, f1], [d0*h, d1*h]);
y = findCubicPolyMinimum(p);

% use only interpolation, no extrapolation permitted
if y < 0 || y > 1 || ~isreal(y) || ~isfinite(y)
    y = 0.5;
end

% safeguards
% if y < 0.001
%     y = 0.001;
% end
% if y>0.9
%     y = 0.9;
% end

% convert from y to t
newt = y*h + t0;
