function x = findCubicPolyMinimum(p)
% findCubicPolyMinimum - find location of the cubic polynomial minimum

% x =  findCubicPolyMinimum(p)  find location of the minimum of cubic
% polynomial whose coefficients are given in 4-vector p, i.e., 
% x = argmin{P(x)}. Where P(x) = p(1)*x^3 + p(2)*x^2 + p(3)*x + p(4)


% convert to a*x^3 + b*x^2 + c*x + d 
% for better readability
a = p(1);
b = p(2);
c = p(3);
d = p(4);


% assuming a /= 0 !!!

if b < 0 % no canellation
    x = (sqrt(b^2-3*a*c) - b)/(3*a);
else % avoid cancellation
    x = -c/(sqrt(b^2 - 3*a*c) + b);
end
    
    
    
    
