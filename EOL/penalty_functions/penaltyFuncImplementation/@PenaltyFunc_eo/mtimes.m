function self = mtimes(self, scalar)

% We need to check which argument is a penalty funciton: it can be either
% the first or second.
if ~isa(self, 'PenaltyFunc_eo')
    [self, scalar] = deal(scalar, self);
end

validateattributes(scalar, {'numeric'}, {'real', 'scalar'});
self.multFactor = scalar * self.multFactor;
