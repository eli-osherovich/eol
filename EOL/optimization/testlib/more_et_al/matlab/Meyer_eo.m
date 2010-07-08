function [fvec, Jac, hessV, hessVfunc, x0_std] = Meyer_eo(x, ~)
%MEYER - Meyer function.

% Number of variables: N = 3
% Number of equations: M = 16
% Standard starting point x0_std = [0.02; 4000; 250] 
% Min function value = 87.9458...
% Argmin x = 


% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [0.02; 4000; 250];

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end

% initialize some constants
t = 45 + 5*(1:16)';
y = [34780, 28610, 23650, 19630, 16370, 13720, 11540, 9744, 8261, ...
     7030, 6005, 5147, 4427, 3820, 3307, 2872]';

fvec = x(1)*exp(x(2)./(t+x(3))) - y;


if nargout > 1 % Jacobian requested
    Jac = zeros(16, 3);
    Jac(:,1) = exp(x(2)./(t+x(3)));
    Jac(:,2) = x(1)*exp(x(2)./(t+x(3)))./(t+x(3));
    Jac(:,3) = -x(1)*x(2)*exp(x(2)./(t+x(3)))./(t+x(3)).^2;
    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:Meyer:HessianNotSupported', ['Hessian or ' ...
                            'Hessian-vector product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
