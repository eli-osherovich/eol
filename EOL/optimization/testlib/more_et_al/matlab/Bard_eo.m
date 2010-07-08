function [fvec, Jac, hessV, hessVfunc, x0_std] = Bard_eo(x, ~)
%BARD - Bard function.

% Number of variables: N = 3
% Number of equations: M = 15
% Standard starting point x0_std = [1; 1; 1] 
% Min function value = 8.21487E-3
% Argmin x = 
% Additional extrema:
% f=17.4286 at x=[0.8406..., -Inf, -Inf]


% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [1; 1; 1];

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end

% initialize some constants
y = [0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 0.58, ...
     0.73, 0.96, 1.34, 2.10, 4.39]';
u = (1:15)';
v = 16 - (1:15)';
w = min(u,v);

fvec = y - (x(1) + u./(v*x(2) + w*x(3)));


if nargout > 1 % Jacobian requested
    Jac = zeros(15, 3);
    Jac(:,1) = -1;
    Jac(:,2) = (v.*u)./(v*x(2) + w*x(3)).^2;
    Jac(:,3) = (w.*u)./(v*x(2) + w*x(3)).^2;
    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:Bard:HessianNotSupported', ['Hessian or ' ...
                            'Hessian-vector product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
