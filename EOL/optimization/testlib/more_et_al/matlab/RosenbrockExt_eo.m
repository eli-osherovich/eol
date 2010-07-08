function [fvec, Jac, hessV, hessVfunc, x0_std] = RosenbrockExt_eo(x, ~, ~, N)
%ROSENBROCKEXT - Extended Rosenbrock function.

% Number of variables: N = any (N must be even) (N=4 by default)
% Number of equations: M = any (M = N)
% Standard starting point x0_std = [-1.2; 1; -1.2; 1; ...; -1.2; 1]
% Min function value = 0
% Argmin x = [1; 1; 1; ...; 1]


% Copyright 2010 Eli Osherovich.

% set N
if 4 > nargin 
    N = 4;
else
    % verify that N is even
    if 0 ~= rem(N,2)
        error('EOL:RosenbrockExt:IllegalArgument','N must be even');
    end
end


% set M
M = N;

% standard starting point
x0_std = repmat([-1.2; 1], N/2, 1);

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end

% initialize constants
I = 1:M/2;

% compute f values
fvec = zeros(M,1);
fvec(2*I-1) = 10*(x(2*I) - x(2*I-1).^2);
fvec(2*I) = 1- x(2*I-1);
if nargout > 1 % Jacobian requested
    Jac = zeros(M,N);
    for i = 1:M/2
        Jac(2*i-1, 2*i) = 10;
        Jac(2*i-1, 2*i-1) = -20*x(2*i-1);
        Jac(2*i, 2*i-1) = -1;
    end
    
    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:RosenbrockExt:HessianNotSupported', ['Hessian or ' ...
                            'Hessian-vector product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
