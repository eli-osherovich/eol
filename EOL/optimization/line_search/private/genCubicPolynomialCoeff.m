function p = genCubicPolynomialCoeff(y, s)
%genCubicPolynomialCoeff - generate the coefficients of the cubic
%polynomial 

%p = genCubicPolynomialCoeff(y, s) - generate a
%vector of legth 4 whose elemets are the coeffiecients in decreasing powers
%of the cubic polyniomial P(x). This polynomial satisfies:
% P(0) = y(0)
% P(1) = y(1)
% P'(0) = s(0)
% P'(1) = s(1)



p(4) = y(1);
p(3) = s(1);
p(2) = 3*(y(2) - y(1)) - 2*s(1) - s(2);
p(1) = 2*(y(1) - y(2))  + s(1) + s(2);
