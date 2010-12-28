classdef LinearOpSum_eo < LinearOp_eo
    properties
        % It should be a list of linear operators, but MATLAB does not
        % allow creation of an array of different linear operators. Hence,
        % the implementation is based on a cell.
        LinOpListCell
    end
    
    
    methods
        function self = LinearOpSum_eo(linopCell)
            % Make sure that all arguments are linops.
            if any(cellfun (@(c) ~isa(c, 'LinearOp_eo'), linopCell))
                error('EOL:LinearOpSum:wrongArgType', ...
                    'Inputs must be of type LinearOp_eo');
            end
            
            % Total number of linops in the chain
            n = numel(linopCell);
            
            % Make sure the linops' deminstions are compatible.
            for i = 2:n
                if linopCell{i}.RangeNumelCurrent ~= linopCell{i-1}.RangeNumelCurrent ||...
                        linopCell{i}.ImageNumelCurrent ~= linopCell{i-1}.ImageNumelCurrent
                    error('EOL:LinearOpSum:dimensionMismatch', ...
                        'PLUS arguments dimensions are not consistent.');
                end
            end
            
            % Create a LinOpSum
            self = self@LinearOp_eo(linopCell{1}.RangeNumelCurrent, linopCell{1}.ImageNumelCurrent);
            self.LinOpListCell = linopCell(:);
        end
        
        function Ax = ApplyForward(self, x)
            Ax = 0;
            for i = numel(self.LinOpListCell):-1:1
                Ax = Ax + self.LinOpListCell{i} * x;
            end
        end
        
        function Ax = ApplyAdjoint(self, x)
            Ax = 0;
            for i = numel(self.LinOpListCell):-1:1
                Ax = Ax + self.LinOpListCell{i}' * x;
            end
        end
    end
end

