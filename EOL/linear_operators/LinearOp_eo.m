classdef LinearOp_eo
    
    
    % Copyright 2010 Eli Osherovich.
    
    properties 
        RangeNumelOriginal
        RangeNumel
        
        ImageNumelOriginal
        ImageNumel
        
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
            self.RangeNumel = rangeNumel;
            
            self.ImageNumelOriginal = imageNumel;
            self.ImageNumel = imageNumel;
        end

        % Toggle the conjugate transpose flag on a call to ctranspose. By
        % defining this function we allow correct treatment of the '
        % (apostrophe) operator.
        function self = ctranspose(self)
            self.AdjointFlag = ~self.AdjointFlag;
            [self.RangeNumel, self.ImageNumel] = deal(self.ImageNumel, self.RangeNumel);
        end
        
        function self = conj(self)
            self.AdjointFlag = ~self.AdjointFlag;
            [self.RangeNumel, self.ImageNumel] = deal(self.ImageNumel, self.RangeNumel);
        end
        
        function self = uminus(self)
            self.MinusFlag = ~ self.MinusFlag;
        end
        
        function newObj = horzcat(varargin)
            newObj = LinearOpMatrix_eo(varargin);
        end
        
        function newObj = vertcat(varargin)
            newObj = LinearOpMatrix_eo(varargin');
        end
        
        function Ax = mtimes(self, x)
            % If x is a linear opertor, the result is a linear operator
            % chain (which is, of course, a linear operator too).
            if isa(x, 'LinearOp_eo')
                Ax = LinearOpChain_eo({self, x});
                return;
            end
            
            % Check that the X's size is consistent with the rangeNumel of
            % the operator.
            if self.RangeNumel ~= numel(x)
                error('EOL:LinearOp:mtimes:wrongDimension', ...
                    'Matrix dimensions must agree.');
            end
            
            % Apply operator in an appropriate mode, either forward or
            % adjoint.
            if self.AdjointFlag
                Ax = ApplyAdjoint(self, x);
            else
                Ax = ApplyForward(self, x);
            end
            
            if self.MinusFlag
                Ax = -Ax;
            end
        end
        
        function newObj = plus(self, x)
            newObj = LinearOpPlus_eo({self, x});
        end
        
    end
    
    methods (Abstract)
        Ax = ApplyAdjoint(self, x)
        Ax = ApplyForward(self, x)
    end
end
