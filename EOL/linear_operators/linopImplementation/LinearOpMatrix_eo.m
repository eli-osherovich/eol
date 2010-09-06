classdef LinearOpMatrix_eo < LinearOp_eo
    
    
    % Copyright 2010 Eli Osherovich.
    
    properties 
        % Matlab does not allow creation of arrays of "different" objects,
        % even if they have a common superclass.
        % Hence, we use a cell here.
        LinopMatrixCell = {}
        AllRangeNumel = []
        AllImageNumel = []
    end
    
    methods

        function self = LinearOpMatrix_eo(linopCell)
            
            % Make sure that all inputs are LinearOp(s).
            if any(cellfun (@(c) ~isa(c, 'LinearOp_eo'), ...
                    linopCell))
                error('EOL:LinearOpMatrix:wrongArgType', ...
                    'Inputs must be of type LinearOp_eo');
            end
            
           
            [m, n] = size(linopCell);
            
            % Run over rows
            allImageNumelTmp = zeros(1, m);
            for i = 1:m
                uniqueVals = unique(cellfun(@(c) c.ImageNumel, linopCell(i,:)));
                if ~all(uniqueVals == uniqueVals(1))
                    error('EOL:LinearOpMatrix:dimensionMismatch', ...
                        'CAT arguments dimensions are not consistent.')
                end
                allImageNumelTmp(i) = uniqueVals(1);
            end
            
            % Run over columns
            allRangeNumelTmp = zeros(1, n);
            for j = 1:size(linopCell, 2)
                uniqueVals = unique(cellfun(@(c) c.RangeNumel, linopCell(:,j)));
                if ~all(uniqueVals == uniqueVals(1))
                    error('EOL:LinearOpMatrix:dimensionMismatch', ...
                        'CAT arguments dimensions are not consistent.')
                end
                allRangeNumelTmp(j) = uniqueVals(1);
            end
            self = self@LinearOp_eo(sum(allRangeNumelTmp), sum(allImageNumelTmp));
            self.LinopMatrixCell = linopCell;
            self.AllRangeNumel = allRangeNumelTmp;
            self.AllImageNumel = allImageNumelTmp;
        end
  
        
        function Ax = ApplyForward(self, x)
            % Split x and Ax into chunks corresponding the rangeNumel and
            % imageNumel of the current linear operators in the array.
            xCell = mat2cell(x(:), self.AllRangeNumel, 1);
            AxCell = cell(numel(self.AllImageNumel), 1);
            for i = 1:size(self.LinopMatrixCell, 1)
                AxRow = 0;
                for j = 1:size(self.LinopMatrixCell, 2)
                    AxRow = AxRow + self.LinopMatrixCell{i,j}*xCell{j};
                end
                AxCell{i} = AxRow;
            end
            Ax = vertcat(AxCell{:});
        end
        
        function Ax = ApplyAdjoint(self, x)
            % Split x and Ax into chunks corresponding the imageNumel and
            % rangeNumel of the current linear operators in the array.
            xCell = mat2cell(x(:), self.AllImageNumel, 1);
            AxCell = cell(numel(self.AllRangeNumel), 1);
            for j = 1:size(self.LinopMatrixCell, 2)
                AxColumn = 0;
                for i = 1:size(self.LinopMatrixCell, 1)
                    AxColumn = AxColumn + self.LinopMatrixCell{i,j}'*xCell{i};
                end
                AxCell{j} = AxColumn;
            end
            Ax = vertcat(AxCell{:});
        end  
        
    end
end
