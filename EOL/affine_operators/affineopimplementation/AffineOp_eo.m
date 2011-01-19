classdef AffineOp_eo < Mapping_eo
    % The class implementing affine tansformation:
    % Ax + b, where A is a linear operator (LinearOp_eo) and b is a vector
    % of appropriate size.
    
    
    % Copyright 2011 Eli Osherovich.
    
    
    properties
        LinOp
        B = 0
    end
    
    methods
        function self = AffineOp_eo(linOp, b)
            
            % Verify that the first argument is a linear operator.
            if ~isa(linOp, 'LinearOp_eo')
                error('EOL:AffineOp:WrongArgType',...
                    'The fist argument must be a linear operator, instead got %s.',...
                    class(linOp));
            end
            
            % Create an object.
            self = self@Mapping_eo(linOp.RangeNumel, linOp.ImageNumel);
            self.LinOp = linOp;
            
            if 2 == nargin
                self.B = b(:);
            end
        end
        
        function Ax = ApplyMapping(self, x)
            Ax = ApplyMapping(self.LinOp, x) + self.B;
        end
        
        function Jx = MultJacobian(self, x, ~)
            Jx = MultJacobian(self.LinOp, x);
        end
        
        function Jcx = MultConjJacobian(self, x, ~)
            Jcx = MultConjJacobian(self.LinOp, x);
        end
    end
end

        
        
        
           
