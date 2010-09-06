classdef LinearOpPlus_eo < LinearOp_eo
    properties
        % It should be a list of linear operators, but MATLAB does not
        % allow creation of an array of different linear operators. Hence,
        % the implementation is based on a cell.
        LinOpListCell
    end
    
    
    methods
        function self = LinearOpPlus_eo(linopCell)
            % Make sure that all arguments are linops.
            if any(cellfun (@(c) ~isa(c, 'LinearOp_eo'), linopCell))
                error('EOL:LinearOpMatrix:wrongArgType', ...
                    'Inputs must be of type LinearOp_eo');
            end
            
            % Total number of linops in the chain
            n = numel(linopCell);
            
            % Make sure the linops' deminstions are compatible.
            for i = 2:n
                if linopCell{i}.RangeNumel ~= linopCell{i-1}.RangeNumel ||...
                        linopCell{i}.ImageNumel ~= linopCell{i-1}.ImageNumel
                    error('EOL:LinearOpPlus:dimensionMismatch', ...
                        'PLUS arguments dimensions are not consistent.');
                end
            end
            
            % Create a LinOpPlus
            self = self@LinearOp_eo(linopCell{1}.RangeNumel, linopCell{1}.ImageNumel);
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

