function [fvec, Jac, hessV, hessVfunc, x0_std] = BrownDennis_eo(x, ...
                                                  ~, M)
%BROWNDENNIS - Brown and Dennis function.

% Number of variables: N = 4
% Number of equations: M = any M>=n (M=20 by default)
% Standard starting point x0_std = [25; 5; -5; 1]
% Min function value = 85822.2 (if M=20)
% Argmin x =


% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [25; 5; -5; 1];

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end

% set number of equations (if not provided)
if 3 > nargin
    M = 20;
end

% initialize some constants
t = (1:M)'/5;

s1 = x(1) + t*x(2) - exp(t);
s2 = x(3) + x(4)*sin(t) - cos(t);

fvec = s1.^2 + s2.^2;


if nargout > 1 % Jacobian requested
    Jac = zeros(M, 4);
    Jac(:, 1) = 2*s1;
    Jac(:, 2) = 2*s1.*t;
    Jac(:, 3) = 2*s2;
    Jac(:, 4) = 2*s2.*sin(t);
    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:BrownDennis:HessianNotSupported', ['Hessian or ' ...
                            'Hessian-vector product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
