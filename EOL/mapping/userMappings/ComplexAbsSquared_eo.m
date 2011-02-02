classdef ComplexAbsSquared_eo < Mapping_eo
    
    methods
        function self = ComplexAbsSquared_eo(rangeNumel)
            if nargin ~= 1
                error('EOL:ComplexAbsSquared_:WrongArgNum', ...
                    'You must provide exactly one argument');
            end
            self = self@Mapping_eo(rangeNumel);
        end
        
        function Ax = ApplyMapping(~, x)
            Ax = abs(x(:)).^2;
        end
        
        function Jx = MultJacobian(~, x, xCurrent)
            Jx = 2*x(:).*conj(xCurrent(:));
        end
        
        function Jcx = MultConjJacobian(~, x, xCurrent)
            Jcx = 2*x(:).*xCurrent(:);
        end
    end
end
