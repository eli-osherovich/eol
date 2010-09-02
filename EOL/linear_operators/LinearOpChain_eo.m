classdef LinearOpChain_eo < LinearOp_eo
    properties
        % It should be a list of linear operators, but MATLAB does not
        % allow creation of an array of different linear operators. Hence,
        % the implementation is based on a cell.
        LinOpListCell
    end
    
    
    methods
        function self = LinearOpChain_eo(linopCell)
            % Make sure that all arguments are linops.
            if any(cellfun (@(c) ~isa(c, 'LinearOp_eo'), linopCell))
                error('EOL:LinearOpMatrix:wrongArgType', ...
                    'Inputs must be of type LinearOp_eo');
            end
            
            % Total number of linops in the chain
            n = numel(linopCell);
            
            % Make sure the linops' deminstions are compatible.
            for i = 2:n
                if linopCell{i}.ImageNumel ~= linopCell{i-1}.RangeNumel
                    error('EOL:LinearOpChain:dimensionMismatch', ...
                        'CAT arguments dimensions are not consistent.');
                end
            end
            
            % Create a LinOpChain
            self = self@LinearOp_eo(linopCell{end}.RangeNumel, linopCell{1}.ImageNumel);
            self.LinOpListCell = linopCell(:);
        end
        
        function Ax = ApplyForward(self, x)
            Ax = x;
            for i = numel(self.LinOpListCell):-1:1
                Ax = self.LinOpListCell{i} * Ax;
            end
        end
        
        function Ax = ApplyAdjoint(self, x)
            Ax = x;
            for i = 1:numel(self.LinOpListCell)
                Ax = self.LinOpListCell{i}' * Ax;
            end
        end
    end
end

