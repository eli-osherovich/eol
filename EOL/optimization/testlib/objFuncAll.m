function [f, ftf, fj, g, x] = objFuncAll(libFile, x, m, n)
%OBJFUNALL - calculate all quantities of a minimization problem.
%
% Usage:
%
%[F, FTF, FJ, G, X] = OBJFUNCALL(LIBFILE, X, M, N)
%
% where
%
% Inputs:
% -------
% LIBFILE - path to a shared library that implements the optimization
% problem (mandatory).
%
% X - current point. Optional argument. If not provided (or empty) the
% default starting point will be used as defined by the optimization
% problem.
% 
% M - number of equations. Optional argument. It should be used only if the
% problem allows this parameter to vary. MATLAB may crash if used
% incorrectly. 
% 
% N - number of variables. Optional argument. It should be used only if the
% problem allows this parameter to vary. MATLAB may crash if used
% incorrectly.
%
% Outputs:
% ---------
% F - a column vector of residuals: F(X).
%
% FTF - squared norm of the residuals vector: F(X)'*F(X).
%
% FJ - Jacobian of F(X).
%
% G - gradient of of FTF(X): 2*FJ'*F
%
% Implementation details
% -----------------------
% Shared library LIBFILE should provide GETFUN_ function. For details see
% Readme file.

% Copyright 2010 Eli Osherovich.


[libDir, libName] = fileparts(libFile);

if ~libisloaded(libName)
    loadlibrary(libFile, fullfile(libDir, 'header.h'));
end


% dummy constants
zero = 0;

% get m and n
[~, nDefault, ~, mDefault, ~, ~, ~, ~, ~] =calllib(libName, 'getfun_',...
    zero, zero, zero, zero, zero, zero, zero, zero, -1) ;

% verify dimensions
if nargin > 2 && ~isempty(m)
    if m ~= mDefault
        warning('EOL:objFunctionAll:DimensionsMismatch', ...
            '%s: Number of equations does not match', libName);
    end
else
    m = mDefault;
end

if nargin > 3 && ~isempty(n)
    if n ~= nDefault
        warning('EOL:objFunctionAll:DimensionsMismatch', ...
            '%s: Number of variables does not match', libName);
    end
else
    n = nDefault;
end


if nargin < 2 || isempty(x)
    % use default x
    x = zeros(n,1);
    [x, ~, ~, ~, ~, ~, ~, ~, ~] = calllib(libName, 'getfun_', x, n, ...
        zero, m, zero, zero, zero, zero, 0) ;
end

if numel(x) ~= n
    error('EOL:objFunctionAll:DimensionsMismatch', ...
        '%s: Number of variables does not match', libName);
end

% preallocate space
f = zeros(m, 1);
ftf = 0;
fj = zeros(m, n);
lfj = m;
g = zeros(n, 1);

% calculate all possible outputs
[~, ~, f, ~, ftf, fj, ~, g, ~] = calllib(libName, 'getfun_', x, n, ...
    f, m, ftf, fj, lfj, g, 1111);
g = 2*g;
