function [fvec, Jac, hessV, hessVfunc, x0_std] = JennrichSampson_eo(x, ...
                                                  ~, M)
%JENNRICHSAMPSON - Jennrich and Sampson function.

% Number of variables: N = 2
% Number of equations: M = any (10 by default)
% Standard starting point x0_std = [0.3; 0.4] 
% Min function value = 124.362 (for M=10)
% Argmin x = [0.2578; 0.2578] (for M=10)



% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [0.3; 0.4];

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end

% compute f values
if 3 > nargin
    M = 10;
end

fvec = zeros(M,1);
for k = 1:M
    fvec(k) = 2 + 2*k - (exp(k*x(1)) + exp(k*x(2)));
end

if nargout > 1 % Jacobian requested
    Jac = zeros(M, 2);
    for k = 1:M
        Jac(k,:) = -k*[exp(k*x(1)), exp(k*x(2))];
    end
    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:JennrichSampson:HessianNotSupported', ['Hessian ' ...
                            'or Hessian-vector product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
