function Ax = applyLinOpAdj(FAstruct, x, idx)


% Copyright 2008-2010 Eli Osherovich.

% Apply linear operator(s) 
if nargin == 3 && isscalar(idx)
    Ax = FAstruct(idx).linop' * x;
else
    Ax = cell(1, length(FAstruct));
    for k = 1:length(Ax)
        Ax{k} = FAstruct(k).linop' * x;
    end
end
