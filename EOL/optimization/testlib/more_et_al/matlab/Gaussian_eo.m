function [fvec, Jac, hessV, hessVfunc, x0_std] = Gaussian_eo(x, ~)
%GAUSSIAN - Gaussian function.

% Number of variables: N = 3
% Number of equations: M = 15
% Standard starting point x0_std = [0.4; 1; 0] 
% Min function value = 1.12793E-8
% Argmin x = 


% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [0.4; 1; 0];

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end

% initialize some constants
t = (8-(1:15)')/2;
y = [0.0009, 0.0044, 0.0175, 0.0540, 0.1295, 0.2420, 0.3521, ...
          0.3989, 0.3521, 0.2420,  0.1295, 0.0540, 0.0175, 0.0044, ...
          0.0009]';

fvec = x(1)*exp(-x(2)*(t - x(3)).^2/2) - y;


if nargout > 1 % Jacobian requested
    Jac = zeros(15, 3);
    Jac(:,1) = exp(-x(2)*(t - x(3)).^2/2);
    Jac(:,2) = -x(1)*exp(-x(2)*(t - x(3)).^2/2).*(t - x(3)).^2/2;
    Jac(:,3) = x(1)*exp(-x(2)*(t - x(3)).^2/2).*(x(2)*(t-x(3)));
    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:Gaussian:HessianNotSupported', ['Hessian or ' ...
                            'Hessian-vector product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
