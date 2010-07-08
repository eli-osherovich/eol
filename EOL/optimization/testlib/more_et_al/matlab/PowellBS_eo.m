function [fvec, Jac, hessV, hessVfunc, x0_std] = PowellBS_eo(x, ~)
%POWELLBS - Powell badly scaled function.

% Number of variables: N = 2
% Number of equations: M = 2
% Standard starting point x0_std = [0; 1] 
% Min function value = 0
% Argmin x = [1.098E-5; 9.106]



% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [0; 1];

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
        1000*x(1)*x(2) - 1
        exp(-x(1)) + exp(-x(2)) - 1.0001
    ];

if nargout > 1 % Jacobian requested
    Jac =...
        [ 
            1000*x(2),     1000*x(1)
            -1/exp(x(1)),   -1/exp(x(2))        
        ];
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:PowellBS:HessianNotSupported', ['Hessian or Hessian-vector ' ...
                            'product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
