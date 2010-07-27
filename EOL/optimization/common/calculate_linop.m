function Ax = calculate_linop(mode, FAstruct, x, idx)


% Copyright 2008-2010 Eli Osherovich.

% do (almost) nothing if x is empty
% if isempty(x)
%     if nargin == 4 && isscalar(idx)
%         Ax = [];
%     else
%         Ax = cell(1, length(FAstruct));
%     end
%     return
% end

% calculte linear operator(s) 
if nargin == 4 && isscalar(idx)
    Ax = FAstruct(idx).linop(mode, x, FAstruct(idx).linop_args{:});
else
    Ax = cell(1, length(FAstruct));
    for k = 1:length(Ax)
        Ax{k} = FAstruct(k).linop(mode, x, FAstruct(k).linop_args{:});
    end
end
