function [fvec, Jac, hessV, hessVfunc, x0_std] = Beal_eo(x, ~)
%BEAL - Beal function.

% Number of variables: N = 2
% Number of equations: M = 3
% Standard starting point x0_std = [1; 1] 
% Min function value = 0
% Argmin x = [3, 0.5]



% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [1; 1];

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
        1.5   - x(1)*(1 - x(2))
        2.25  - x(1)*(1 - x(2)^2)
        2.625 - x(1)*(1 - x(2)^3)
    ];

if nargout > 1 % Jacobian requested
    Jac =...
        [ 
            -1 + x(2),     x(1)
            -1 + x(2)^2,   2*x(1)*x(2)
            -1 + x(2)^3,   3*x(1)*x(2)^2
        ];
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:Beale:HessianNotSupported', ['Hessian or Hessian-vector ' ...
                            'product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
