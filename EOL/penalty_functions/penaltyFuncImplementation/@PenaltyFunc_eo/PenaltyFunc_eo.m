classdef PenaltyFunc_eo
    
    
    
    % Copyright 2010 Eli Osherovich.
    
    
    
    properties (Access = protected)
        multFactor = 1
    end
    
    methods (Abstract)
        [val, grad, hessMultVecorFunc] = doCalculations(self, x)
    end
end
