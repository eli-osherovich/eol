function [fvec, Jac, hessV, hessVfunc, x0_std] = BrownBS_eo(x, ~)
%BROWNBS - Brown badly scaled function.

% Number of variables: N = 2
% Number of equations: M = 3
% Standard starting point x0_std = [1; 1] 
% Min function value = 0
% Argmin x = [1E6, 2E-6]



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
        x(1) - 10^6
        x(2) - 2*10^-6
        x(1)*x(2) - 2
    ];

if nargout > 1 % Jacobian requested
    Jac =...
        [ 
            1,      0
            0,      1
            x(2),   x(1)
        ];
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:BrownBS:HessianNotSupported', ['Hessian or Hessian-vector ' ...
                            'product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
