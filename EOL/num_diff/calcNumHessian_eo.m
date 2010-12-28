function NHess = calcNumHessian_eo(x, func, mode)
% CALCNUMHESSIAN - calculate Hessian numerically.
% NHESS = CALCNUMHESSIAN(X, FUNC, MODE) - calculate numerical
% approximation to the Hessian of a function FUNC at point X.
% 
% FUNC must be a function (handle) that accepts a single input X
% and returns two outputs:
% [FVAL, GRAD] = FUNC(X) with FVAL being the function value and
% GRAD the function gradient at X.
%
% The optional input MODE can be either 'forward' for forward difference,
% 'backward' for backward difference, 'central' for central difference
% (default), or 'precise' for central difference with Richardson's
% extrapolation (may be slow).



% Copyright 2007-2010 Eli Osherovich.



% Number of variables.
N = numel(x);

% Preallocate space for the Hessian.
NHess = zeros(N);

% Use central difference by default.
if nargin < 3
    mode = 'central';
end

% Calculate partial derivatives df/dx_i
for i = 1:N
    NHess(:,i) = calcNumJacobian_eo(x, @(x) funcWrapper(func, x), mode);
end


function grad = funcWrapper(x, func)
% FUNCWRAPPER - a function wrapper that makes GRAD to be the output
% of FUNC.
[~, grad] = func(x);
