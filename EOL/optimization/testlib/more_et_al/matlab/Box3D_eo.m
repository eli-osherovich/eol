function [fvec, Jac, hessV, hessVfunc, x0_std] = Box3D_eo(x, ~, M)
%BOX3D - Box three-dimensional function.

% Number of variables: N = 3
% Number of equations: M = any M>=N (M=10 by default)
% Standard starting point x0_std = [0; 10; 20] 
% Min function value = 0
% Argmin x = [1; 10; 1], [10; 1; -1] and whenever x(2)=x(1) and x(3)=0


% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [0; 10; 20];

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end

% number of equations
if 3 > nargin
    M = 10;
end

% initialize some constants
t = 0.1*(1:M)';

fvec =  exp(-t*x(1)) - exp(-t*x(2)) - x(3)*(exp(-t) - exp(-10*t));


if nargout > 1 % Jacobian requested
    Jac = zeros(M, 3);
    Jac(:,1) = -t.*exp(-t*x(1));
    Jac(:,2) = t.*exp(-t*x(2));
    Jac(:,3) = -exp(-t) + exp(-10*t);
    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:Box3D:HessianNotSupported', ['Hessian or ' ...
                            'Hessian-vector product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
