function [fvec, Jac, hessV, hessVfunc, x0_std] = PowellSingular_eo(x, ~)
%POWELLSINGULAR - Powell singular function.

% Number of variables: N = 4
% Number of equations: M = 4
% Standard starting point x0_std = [3; -1; 0; 1]
% Min function value = 0
% Argmin x = [0; 0; 0; 0]


% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [3; -1; 0; 1];

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
        x(1) + 10*x(2)
        sqrt(5)*(x(3) - x(4))
        (x(2) - 2*x(3))^2
        sqrt(10)*(x(1) - x(4))^2
    ];


if nargout > 1 % Jacobian requested
    Jac =...
        [
            1,                      10,                0,          0
            0,                      0,                 sqrt(5),    -sqrt(5)
            0,                      2*(x(2) - 2*x(3)), -4*(x(2)-2*x(3)), 0
            2*sqrt(10)*(x(1)-x(4)), 0,                 0,       -2*sqrt(10)*(x(1)-x(4))
        ];
    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:PowellSingular:HessianNotSupported', ['Hessian or ' ...
                            'Hessian-vector product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
