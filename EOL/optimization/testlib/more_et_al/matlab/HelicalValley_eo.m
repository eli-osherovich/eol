function [fvec, Jac, hessV, hessVfunc, x0_std] = HelicalValley_eo(x, ~)
%HELICALVALLEY - Helical valley function.

% Number of variables: N = 3
% Number of equations: M = 3
% Standard starting point x0_std = [-1; 0; 0] 
% Min function value = 0
% Argmin x = [1; 0; 0]


% Implementation details
% ----------------------
% It seems that x(2) is assumed to be positive.
% If this assumption is not correct the function theta(x1,x2)
% should be fixed in order to remain continuous.


% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [-1; 0; 0];

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end


fvec = ...
    [
        10*(x(3) - 10*theta(x(1), x(2)))
        10*(sqrt(x(1)^2 + x(2)^2) - 1)
        x(3)
    ];

if nargout > 1 % Jacobian requested
    Jac = zeros(3);
    
    Jac(1,1) = 50*x(2)/(pi*(x(1)^2 + x(2)^2));
    Jac(1,2) = -50*x(1)/(pi*(x(1)^2 + x(2)^2));
    Jac(1,3) = 10;
    
    Jac(2,1) = 10*x(1)/(x(1)^2 + x(2)^2)^(1/2);
    Jac(2,2) = 10*x(2)/(x(1)^2 + x(2)^2)^(1/2); 

    Jac(3,3) = 1;

    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:HelicalValley:HessianNotSupported', ['Hessian ' ...
                            'or Hessian-vector product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end

function th = theta(x1, x2)
if x1 >= 0 
    th = atan(x2/x1)/(2*pi);
else
    th = atan(x2/x1)/(2*pi) + 0.5;
end
