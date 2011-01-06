classdef Matrix_eo < LinearOp_eo
    
    properties
        Matrix
    end
    
    methods 
        function self = Matrix_eo(matrix)
            
            % Make sure we got a matrix.
            validateattributes(matrix, {'numeric'}, {'2d'});
            
            [m,n] = size(matrix);
            
            self = self@LinearOp_eo(n, m);
            self.Matrix = matrix;
        end
        
        function Ax = ApplyForward(self, x)
            Ax = self.Matrix * x;
        end
        
        function Ax = ApplyAdjoint(self, x)
            Ax = self.Matrix' * x;
        end
    end
end
