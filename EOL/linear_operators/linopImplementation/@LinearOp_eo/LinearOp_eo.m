classdef LinearOp_eo < Mapping_eo
    
    
    % Copyright 2010 Eli Osherovich.
    
    properties 
        % The *Original variables keep the range and image dimensionality
        % at the creation time. Applying conjugate (adjoint) operations
        % interchange the above dimensionalities.
        
        AdjointFlag = false
    end
  
    methods
        function self = LinearOp_eo(rangeNumel, imageNumel)
            if nargin == 1
                % If imageNumel is not specified assume it is equal to
                % rangeNumel.
                imageNumel = rangeNumel;
            end
                
            self = self@Mapping_eo(rangeNumel, imageNumel);
        end
    end
    
    methods 
        function Ax = ApplyMapping(self, x)
            if self.AdjointFlag
                Ax = ApplyAdjoint(self, x);
            else
                Ax = ApplyForward(self, x);
            end
        end
        
        function Jx = MultJacobian(self, x, ~)
            if self.AdjointFlag
                Jx = ApplyAdjoint(self, x);
            else
                Jx = ApplyForward(self, x);
            end
        end
        
        function Jcx = MultConjJacobian(self, x, ~)
            if self.AdjointFlag
                Jcx = ApplyForward(self, x);
            else
                Jcx = ApplyAdjoint(self, x);
            end
        end
        
    end
    
    methods (Abstract)
        Ax = ApplyAdjoint(self, x)
        Ax = ApplyForward(self, x)
    end
end
