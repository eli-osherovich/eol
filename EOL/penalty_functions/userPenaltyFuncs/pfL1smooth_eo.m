classdef pfL1smooth_eo < PenaltyFunc_eo
    % Penalize (weighted) smooth L1 norm approximation of x-x0: 
    % w.*sqrt((x-x0)^2 + e)
    
    
    % Copyright 2010 Eli Osherovich.
    
    
    properties
        e = 0.001
        x0 = 0
        w = 1
    end

    methods
        function self = pfL1_eo(e, x0, w)
            if nargin > 0
                validateattributes(e, {'numeric'}, {'real', 'positive'});
                self.e = e;
                
                if nargin > 1
                    validateattributes(x0, {'numeric'}, {'finite'});
                    self.x0 = x0(:);
                    if nargin > 2
                        validateattributes(w, {'numeric'}, {'real', 'nonnegative'});
                        self.w = w(:);
                    end
                end
            end
        end
        
        function [val, grad, hessMultVecorFunc] = doCalculations(self, x)
            tmp = sqrt((x - self.x0).^2 + self.e);
            val = sum(self.w .* (tmp - self.e));
            
            if nargout > 1, % gradient requested
                grad = self.w .* (x-self.x0)./tmp;
                if nargout > 2, % Hessian mult. function requested
                    hessMultVecorFunc = @hessMult;
                end
            end
            
            function hessV = hessMult(v)
                hessV =  self.e .*v ./ tmp.^3;
            end
        end
    end
end
