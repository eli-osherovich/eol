classdef pfL1_eo < PenaltyFunc_eo
    % Penalize the weighted sum: sum(W.*abs(X-x0))
    
    
    % Copyright 2010 Eli Osherovich.
    
    
    properties
        x0 = 0
        w = 1
    end

    methods
        function self = pfSum_eo(x0, w)
            switch nargin
                case 0
                    % Do nothing
                case 1
                    validateattributes(x0, {'numeric'}, {'finite'});
                    self.x0 = x0(:);
                case 2
                    validateattributes(x0, {'numeric'}, {'finite'});
                    validateattributes(w, {'numeric'}, {'real', 'nonnegative'});
                    self.x0 = x0(:);
                    self.w = w(:);
            end
        end
        
        function [val, grad, hessMultVecorFunc] = doCalculations(self, x)
            val = sum(self.w .* abs(x-self.x0));
            
            if nargout > 1, % gradient requested
                grad = self.w .* sign(x);
                if nargout > 2, % Hessian mult. function requested
                    hessMultVecorFunc = @hessMult;
                end
            end
            
            function hessV = hessMult(v)
                hessV =  0 * v;
            end
        end
    end
end
