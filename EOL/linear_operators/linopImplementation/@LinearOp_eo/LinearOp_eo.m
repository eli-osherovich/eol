classdef LinearOp_eo
    
    
    % Copyright 2010 Eli Osherovich.
    
    properties 
        RangeNumelOriginal
        RangeNumelCurrent
        
        ImageNumelOriginal
        ImageNumelCurrent
        
        AdjointFlag = false
        MinusFlag = false
    end
  
    methods
        function self = LinearOp_eo(rangeNumel, imageNumel)
            if nargin == 1
                % If imageNumel is not specified assume it is equal to
                % rangeNumel.
                imageNumel = rangeNumel;
            end
                
            validateattributes(rangeNumel, {'numeric'}, {'integer', ...
                'positive', 'scalar'});
            validateattributes(imageNumel, {'numeric'}, {'integer', ...
                'positive', 'scalar'});
            
            self.RangeNumelOriginal = rangeNumel;
            self.RangeNumelCurrent = rangeNumel;
            
            self.ImageNumelOriginal = imageNumel;
            self.ImageNumelCurrent = imageNumel;
        end
    end
    
    methods (Abstract)
        Ax = ApplyAdjoint(self, x)
        Ax = ApplyForward(self, x)
    end
end
