classdef PenaltyFunc_eo
    
    
    
    % Copyright 2010 Eli Osherovich.
    
    
    
    properties (Access = private)
        multFactor = 1
        minusFlag = false
    end
    
    methods (Abstract)
        [val, grad, hessMultVecorFunc] = doCalculations(self, x)
    end
end
