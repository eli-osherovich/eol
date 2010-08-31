classdef LinearOp_eo
    
    properties (Access = protected)
        RangeSize
        ImageSize
        AdjointFlag = false
    end
    
    methods        
        function self = LinearOp_eo(rangeSize, imageSize)
            if nargin < 1
                error('EOL:LinearOp:WrongInput', ['Not enough input ' ...
                    'arguments']);
            end
            
            
            % Inputs must be integer positive vectors.
            validateattributes(rangeSize, {'numeric'}, {'integer', ...
                'positive', 'vector'});
            self.RangeSize = rangeSize;
            
            % Set image shape (size).
            if nargin > 1
                self.ImageSize = imageSize;
            else
                self.ImageSize = rangeSize;
            end
            
        end
        
        % Toggle the conjugate transpose flag on a call to ctranspose. By
        % defining this function we allow correct treatment of the '
        % (apostrophe) operator.
        function self = ctranspose(self)
            self.AdjointFlag = ~self.AdjointFlag;
        end
    end
    
    methods (Abstract)
        Ax = mtimes(self, x)
    end
end
