function valid = isValid(x, allowComplex)
% ISVALID - check vector validity.


% Copyright 2010 Eli Osherovich.

% Do not allow complex values by default.
if nargin < 2 || isempty(allowComplex)
    allowComplex = false;
end

valid = (allowComplex || isreal(x)) && all(isfinite(x(:)));
