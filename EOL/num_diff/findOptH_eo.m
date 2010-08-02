function h = findOptH_eo(x, mode)
% FINDOPTH - find optimal h for numerical differentiation. 
% 
% H = FINDOPTH(X) - Find optimal H for central difference approximation.
% H = FINDOPTH(X, MODE) - Find optimal H for central/forward/backward
% difference approximation.
%
%
% Inputs:
% --------------------
% X - a real numeric array.
% MODE - a string indicating the approximation mode (can be 'forward',
% 'backward', or 'central').
% 
%
% Outputs:
% --------------------
% H -  a double array whose shape matches that of X.



% Copyright 2010 Eli Osherovich.


 
% Make sure x is real numeric.
validateattributes(x, {'numeric'}, {'real'});

% Central difference is used by default.
if nargin < 2
    mode = 'central';
end

switch lower(mode)
  case 'forward'
    h = eps^(1/2);
    % Update h such that x+h is exactly representable.
    h = (x + h) - x;
  case 'backward'
    h = eps^(1/2);
    % Update h such that x-h is exactly representable.
    h = x - (x - h);
  case {'central', 'precise'}
    h = eps^(1/3);
    % Update h such that x+-h is exactly representable.
    h = (abs(x) + h) - abs(x);
  otherwise
    error('EOL:FINDOPTH:WrongInput', ['Unknown mode: %s it must be either ' ...
                        'forward, backward, central, or precise'], mode);
end

% Make sure that h is large enough such that x+h ~= x
h = max(h, eps(x));
