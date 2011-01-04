function Ax = applyMapping(FAstruct, x, idx)


% Copyright 2008-2011 Eli Osherovich.

% Apply mapping(s).
if nargin == 3 && isscalar(idx)
    Ax = ApplyMapping(FAstruct(idx).mapping, x);
else
    Ax = cell(1, length(FAstruct));
    for k = 1:length(Ax)
        Ax{k} = ApplyMapping(FAstruct(k).mapping, x);
    end
end
