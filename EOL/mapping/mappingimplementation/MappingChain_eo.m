classdef MappingChain_eo < Mapping_eo
    
    
    
    % Copyright 2010-2011 Eli Osherovich.
 
    
     
    properties
        % A list of mapping is kept in a cell array, since MATLAB does not
        % allow array of different (subclasses of) Mapping(s).
        % Mappings are stored in their lexical order, i.e., given that the
        % list equals {M1 M2 M3} the resulting mapping would be M1*M2*M3*x.
        MappingListCell
    end
    
    methods
        function self = MappingChain_eo(mappingCell)
            
            % Total number of mappings in the chain.
            n = numel(mappingCell);
            
            % Make sure the list is not empty.
            if 0 == n
                error('EOL:MappingChain:WrongArgNum', ...
                    'You should provide at least one mapping');
            end
            
            % Make sure that all arguments are Mappings(s).
            if any(cellfun (@(c) ~isa(c, 'Mapping_eo'), mappingCell))
                error('EOL:MappingChain_eo:wrongArgType', ...
                    'Inputs must be of type Mapping_eo');
            end
            
            % Make sure that the dimensions are compatible.
            for i = 2:n
                if mappingCell{i}.ImageNumel ~= ...
                        mappingCell{i-1}.RangeNumel
                    error('EOL:MappingChain:dimensionMismatch', ...
                        'Dimensions are inconsistent between mappings %d and %d', i, i-1);
                end
            end
            
            % Create MappingChain.
            self = self@Mapping_eo(mappingCell{end}.RangeNumel, ...
                mappingCell{1}.ImageNumel);
            
            % We could use a simple construction like the following
            % self.MappingListCell = mappingCell(:);
            % However, this is very inefficient: one creates many
            % MappingChanin(s) for an expression like M1*M2*...*M10.
            % Furthemore, evaluating such a chain will result in many
            % *nested* calls which are not necessary.
            tmpCell = {};
            for i = 1:n
                if isa(mappingCell{i}, 'MappingChain_eo')
                    tmpCell = [tmpCell, mappingCell{i}.MappingListCell];
                else
                    tmpCell = [tmpCell mappingCell(i)];
                end
            end
            self.MappingListCell = tmpCell;
            
        end
        
        function Ax = ApplyMapping(self, x)
            Ax = x;
            for i = numel(self.MappingListCell):-1:1
                Ax = ApplyMapping(self.MappingListCell{i}, Ax);
            end
        end
        
        function Jx = MultJacobian(self, x, xCurrent)
            Jx = x;
            for i = numel(self.MappingListCell):-1:1
                Jx = MultJacobian(self.MappingListCell{i}, ...
                    Jx, xCurrent);
            end
        end
            
        function Jcx = MultConjJacobian(self, x, xCurrent)
            Jcx = x;
            for i = 1:numel(self.MappingListCell)
                Jcx = MultConjJacobian(self.MappingListCell{i}, ...
                    Jcx, xCurrent);
            end
        end
    end
end
