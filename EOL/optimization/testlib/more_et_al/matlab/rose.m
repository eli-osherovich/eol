function [fvec, Jac, hessV, hessVfunc, x0_std] = rose(x, ~)
%ROSE - Rosenbrock function.

% Number of equations = 2
% Number of variables = 2
% Min function value = 0
% Argmin x = [1;1]
% Standard starting point x0_std = [-1.2; 1] 


% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [-1.2; 1];

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end

% compute f values
fvec = ...
    [
        10 * (x(2) - x(1)^2);
        1 - x(1)
    ];

if nargout > 1 % Jacobian requested
    Jac =...
        [ 
            -20*x(1),  10
            -1,        0 
        ];
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:Rose:HessianNotSupported', ['Hessian or Hessian-vector ' ...
                            'product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
