function self = mtimes(self, scalar)

% We need to check which is argument is a penalty funciton.
if ~isa(self, 'PenaltyFunc_eo')
    [self, scalar] = deal(scalar, self);
end

validateattributes(scalar, {'numeric'}, {'real', 'scalar'});
self.multFactor = scalar;
