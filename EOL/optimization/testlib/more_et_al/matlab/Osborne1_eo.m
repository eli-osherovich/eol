function [fvec, Jac, hessV, hessVfunc, x0_std] = Osborne1_eo(x, ~)
%OSBORNE1 - Osborne 1 function.

% Number of variables: N = 5
% Number of equations: M = 33;
% Standard starting point x0_std = [0.5; 1.5; -1; 0.01; 0.02]
% Min function value = 5.46489E-5
% Argmin x =


% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [0.5; 1.5; -1; 0.01; 0.02];

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end

% initialize some constants
t = 10*((1:33)' - 1);
y = [0.844; 0.908; 0.932; 0.936; 0.925; 0.908; 0.881; 0.850; 0.818; ...
     0.784; 0.751; 0.718; 0.685; 0.658; 0.628; 0.603; 0.580; 0.558; ...
     0.538; 0.522; 0.506; 0.490; 0.478; 0.467; 0.457; 0.448; 0.438; ...
     0.431; 0.424; 0.420; 0.414; 0.411; 0.406];


fvec = y - (x(1) + x(2)*exp(-t*x(4)) + x(3)*exp(-t*x(5)));


if nargout > 1 % Jacobian requested
    Jac = zeros(33, 5);
    Jac(:, 1) = -1; 
    Jac(:, 2) = -exp(-t*x(4));
    Jac(:, 3) = -exp(-t*x(5));
    Jac(:, 4) = x(2)*t.*exp(-t*x(4));
    Jac(:, 5) = x(3)*t.*exp(-t*x(5));
    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:Osborne1:HessianNotSupported', ['Hessian or ' ...
                            'Hessian-vector product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
