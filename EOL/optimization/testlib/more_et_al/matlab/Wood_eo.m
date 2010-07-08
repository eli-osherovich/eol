function [fvec, Jac, hessV, hessVfunc, x0_std] = Wood_eo(x, ~)
%WOOD - Wood function.

% Number of variables: N = 4
% Number of equations: M = 6
% Standard starting point x0_std = [-3; -1; -3; -1]
% Min function value = 0
% Argmin x = [1; 1; 1; 1]


% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [-3; -1; -3; -1];

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end

fvec =...
    [
        10*(x(2)-x(1)^2)
        1 - x(1)
        sqrt(90)*(x(4) - x(3)^2)
        1 - x(3)
        sqrt(10)*(x(2) + x(4) - 2)
        sqrt(10)*(x(2) - x(4))
    ];


if nargout > 1 % Jacobian requested
    Jac =...
        [
           -20*x(1),   10,       0,               0
           -1,         0,        0,               0
            0,         0,       -2*sqrt(90)*x(3), sqrt(90)
            0,         0,       -1,               0
            0,         sqrt(10), 0,               sqrt(10)
            0,         sqrt(10), 0,               -sqrt(10)
        ];
    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:Wood:HessianNotSupported', ['Hessian or ' ...
                            'Hessian-vector product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
