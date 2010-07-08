function [fvec, Jac, hessV, hessVfunc, x0_std] = Gulf_eo(x, ~, M)
%GULF- Gulf research and development function.

% Number of variables: N = 3
% Number of equations: M = any in the interval [n,100] (100 by default)
% Standard starting point x0_std = [5; 2.5; 0.15] 
% Min function value = 0
% Argmin x = [50; 25; 1.5]


% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [5; 2.5; 0.15];

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end

% number of equations
if 3 > nargin
    M = 100;
end

% initialize some constants
t = (1:M)'/100;
y = 25 + (-50*log(t)).^(2/3);
z = M*y.*(1:M)'; 

tmp = exp(-(abs(z*x(2)).^x(3))/x(1));

fvec =  tmp - t;


if nargout > 1 % Jacobian requested
    Jac = zeros(M, 3);
    Jac(:,1) = tmp .* abs(z*x(2)).^x(3)/x(1)^2;
    Jac(:,2) = -tmp .* abs(z*x(2)).^x(3)*x(3)/(x(1)*x(2));
    Jac(:,3) = -tmp .* abs(z*x(2)).^x(3).*log(abs(z*x(2)))/x(1);
    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:Gulf:HessianNotSupported', ['Hessian or ' ...
                            'Hessian-vector product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
