classdef MappingMatrix_eo < Mapping_eo
    
    
    
    % Copyright 2010-2012 Eli Osherovich.
    
    
    
    properties 
        % Matlab does not allow creation of arrays of "different" objects,
        % even if they have a common superclass.
        % Hence, we use a cell here.
        MappingMatrixCell = {}
        AllRangeNumel = []
        AllImageNumel = []
    end
    
    methods

        function self = MappingMatrix_eo(mappingCell)
            
            % Make sure that all inputs are Mapping(s).
            if any(cellfun (@(c) ~isa(c, 'Mapping_eo'), ...
                    mappingCell))
                error('EOL:MappingMatrix:wrongArgType', ...
                    'Inputs must be of type Mapping_eo');
            end
            
            
            % The size of the new matrix (of mappings).
            [m, n] = size(mappingCell);
            
            % Run over all rows and check that all mappings in the current
            % row have the same number of elements in their image.
            allImageNumelTmp = zeros(1, m);
            for i = 1:m
                uniqueVals = unique(cellfun(@(c) c.ImageNumel, mappingCell(i,:)));
                if ~all(uniqueVals == uniqueVals(1))
                    error('EOL:MappingMatrix:dimensionMismatch', ...
                        'CAT arguments dimensions are not consistent.')
                end
                allImageNumelTmp(i) = uniqueVals(1);
            end
            
            % Run over all columns and check that all mappings in the
            % current column have the same number of elements in their
            % range.
            allRangeNumelTmp = zeros(1, n);
            for j = 1:size(mappingCell, 2)
                uniqueVals = unique(cellfun(@(c) c.RangeNumel, mappingCell(:,j)));
                if ~all(uniqueVals == uniqueVals(1))
                    error('EOL:MappingMatrix:dimensionMismatch', ...
                        'CAT arguments dimensions are not consistent.')
                end
                allRangeNumelTmp(j) = uniqueVals(1);
            end
            self = self@Mapping_eo(sum(allRangeNumelTmp), sum(allImageNumelTmp));
            self.MappingMatrixCell = mappingCell;
            self.AllRangeNumel = allRangeNumelTmp;
            self.AllImageNumel = allImageNumelTmp;
        end
  
        
        function Ax = ApplyMapping(self, x)
            % Split x and Ax into chunks corresponding the rangeNumel and
            % imageNumel of the current mapping in the array.
            xCell = mat2cell(x(:), self.AllRangeNumel, 1);
            AxCell = cell(numel(self.AllImageNumel), 1);
            
            [m, n] = size(self.MappingMatrixCell);
            
            for i = 1:m
                AxRow = 0;
                
                for j = 1:n
                    AxRow = AxRow + self.MappingMatrixCell{i,j} * xCell{j};
                end
                
                AxCell{i} = AxRow;
            end
            
            Ax = vertcat(AxCell{:});
        end
        
        
        
        function Jx = MultJacobian(self, x, xCurrent)
            % Split x and Jx into chunks corresponding the rangeNumel and
            % imageNumel of the mappings in the array.
            xCell = mat2cell(x(:), self.AllRangeNumel, 1);
            JxCell = cell(numel(self.AllImageNumel), 1);
            
            [m, n] = size(self.MappingMatrixCell);
            
            for i = 1:m
                JxRow = 0;
                
                for j = 1:n
                    JxRow = JxRow + ...
                        MultiplyJacobian(self.MappingMatrixCell{i,j},...
                        xCell{j}, xCurrent);
                end
                JxCell{i} = JxRow;
            end
            Jx = vertcat(JxCell{:});
        end

        
        function Jcx = MultConjJacobian(self, x, xCurrent)
            % Split x and Jcx into chunks corresponding the imageNumel and
            % rangeNumel of the mappings in the array.
            xCell = mat2cell(x(:), self.AllImageNumel, 1);
            JcxCell = cell(numel(self.AllRangeNumel), 1);
            
            [m, n] = size(self.MappingMatrixCell);
            
            for j = 1:n
                JcxColumn = 0;
                
                for i = 1:m
                    JcxColumn = JcxColumn + ...
                        MultConjJacobian(self.MappingMatrixCell{i,j}, ...
                        xCell{i}, xCurrent);
                end
                JcxCell{j} = JcxColumn;
            end
            Jcx = vertcat(JcxCell{:});
        end  
        
        
    end
end
