classdef pfDistanceL2Squared_eo < PenaltyFunc_eo
    % Penalize the weighted squared distance  W.*||X - X0||^2
    
    
    % Copyright 2010 Eli Osherovich.
    
    
    properties
        x0 = 0
        w = 1
    end

    methods
        function self = pfDistanceL2Squared_eo(x0, w)
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
            val = sum(self.w .* abs(x-self.x0).^2);
            
            if nargout > 1, % gradient requested
                grad = 2 * self.w .* (x - self.x0);
                if nargout > 2, % Hessian mult. function requested
                    hessMultVecorFunc = @hessMult;
                end
            end
            
            function hessV = hessMult(v)
                hessV = 2 * self.w .* v;
            end
        end
    end
end
