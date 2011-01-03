function Ax = mtimes(self, x)
    
    
    
    % Copyright 2011 Eli Osherovich.
    
    
    % If x is a mapping, the result is a chain of mappings.
    if isa(x, 'Mapping_eo')
        Ax = MappingChain_eo({self, x});
    elseif isnumeric(x)
        Ax = ApplyMapping(self, x);
    elseif iscell(x)
        % Preallocate space.
        Ax = cell(size(x));
        for i = 1:numel(x)
            % Skip if empty.
            if isempty(x{i})
                continue;
            end
            Ax{i} = ApplyMapping(self, x);
        end
    else
        error('EOL:Mapping:mtimes:WrongArgType', ...
            ['x must be a mapping, numeric vector, or cell array.',...
            '\nInstead got %s.'], class(x));
    end
