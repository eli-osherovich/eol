classdef DiagMatrix_eo < LinearOp_eo
    
    properties
        Diagonal
    end
    
    methods 
        function self = DiagMatrix_eo(diagonal)
            
            % Make sure we got a matrix.
            validateattributes(diagonal, {'numeric'}, {'vector'});
            
            n = numel(diagonal);
            
            self = self@LinearOp_eo(n, n);
            self.Diagonal = diagonal;
        end
        
        function Ax = ApplyForward(self, x)
            Ax = self.Diagonal .* x;
        end
        
        function Ax = ApplyAdjoint(self, x)
            Ax = conj(self.Diagonal) .* x;
        end
    end
end
