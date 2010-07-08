function [fvec, Jac, hessV, hessVfunc, x0_std] = KowalikOsborne_eo(x, ~)
%KOWALIKOSBORN - Kowalik and Osborne function.

% Number of variables: N = 4
% Number of equations: M = 11
% Standard starting point x0_std = [0.25; 0.39; 0.415; 0.39]
% Min function value = 3.07505E-4
% Argmin x = 
% Additional extrema:
% f=1.02734 at x=[Inf; -14.07..., -Inf, -Inf]


% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [0.25; 0.39; 0.415; 0.39];

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end


% initialize some constants
y = [0.1957; 0.1947; 0.1735; 0.1600; 0.084; 0.0627; 0.0456; 0.0342; ...
     0.0323; 0.0235; 0.0246;];
u = [4.0000; 2.0000; 1.000; 0.5000; 0.2500; 0.1670; 0.1250; 0.1000; ...
     0.083; 0.0714; 0.0625];

num =  x(1)*(u.^2 + u*x(2));
denom = (u.^2 + u*x(3) + x(4));
fvec = y - num./denom;


if nargout > 1 % Jacobian requested
    Jac = zeros(11, 4);
    Jac(:, 1) = -(u.^2 + u*x(2))./denom;
    Jac(:, 2) = -x(1)*u./denom;
    Jac(:, 3) = num.*u./denom.^2;
    Jac(:, 4) = num./denom.^2;
    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:KowalikOsborne:HessianNotSupported', ['Hessian or ' ...
                            'Hessian-vector product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
