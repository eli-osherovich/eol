classdef pfGEQL1_eo < PenaltyFunc_eo
    
    
    
    % Copyright 2011  Eli Osherovich.
    
    
    properties (Access = private)
        threshold = 0
        w = 1
    end
    
    methods
        function self = pfGEQL1_eo(thr, w)
            switch nargin
                case 0
                    % Do nothing.
                    % This is the default case: x >= 0.
                    
                case 1
                    validateattributes(thr, {'numeric'}, {'real', 'finite'})
                    self.threshold = thr(:);
                case 2
                    validateattributes(thr, {'numeric'}, {'real', 'finite'});
                    validateattributes(w, {'numeric'}, {'real', 'nonnegative'});
                    self.threshold = thr(:);
                    self.w = w(:);
            end
        end
            
        function [val, grad, hessMultVecorFunc] = doCalculations(self, x)
           
            % Find violations.
            difference = x - self.threshold;
            violIdx = difference < 0;
            
            
            % Parameters' values at violating indices. 
            if isscalar(self.w)
                wV = self.w;
            else
                wV = self.w(violIdx);
            end
            differenceV = difference(violIdx);
            
            val = -sum( wV .* differenceV);
            
            if nargout > 1 % gradient requested
                grad = zeros(size(x));
                grad(violIdx) = -wV;
                
                if nargout > 2 % Hessian mult. function requested
                    hessMultVecorFunc = @hessMult;
                end
            end
            
            function hessV = hessMult(v)
                hessV = 0*v;
            end
                
                    
        end
            
    end
end
