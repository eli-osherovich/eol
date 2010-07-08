function [fvec, Jac, hessV, hessVfunc, x0_std] = FreudensteinRoth_eo(x, ~)
%FREUDENSTEINROTH - Freudenstein and Roth function.

% Number of variables: N = 2
% Number of equations: M = 2
% Standard starting point x0_std = [0.5; -2]
% Min function value = 0
% Argmin x = [5; 4]
% Additional extrema:
% f=48.9842 at x=[11.41...; -0.8968...]


% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [0.5; -2];

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
        -13 + x(1) + ((5 - x(2))*x(2) - 2)*x(2);
        -29 + x(1) + ((x(2) + 1)*x(2) - 14)*x(2);
    ];

if nargout > 1 % Jacobian requested
    Jac =...
        [
            1, -3*x(2)^2 + 10*x(2) - 2
            1, 3*x(2)^2 + 2*x(2) - 14
        ];
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:FreudensteinRoth:HessianNotSupported', ['Hessian ' ...
                            'or Hessian-vector product are not supported.']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
