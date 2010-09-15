classdef pfGEQL2Squared_eo < PenaltyFunc_eo
    
    
    
    % Copyright 2008-2010  Eli Osherovich.
    
    
    properties (Access = private)
        threshold = 0
        w = 1
    end
    
    methods
        function self = pfProjGEQL2Squared_eo(thr, w)
            switch nargin
                case 0
                    % Do nothing
                    % This is the default case: the real part of x
                    % projected onto the non-negative orthant.
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
            
            val = sum( wV .* differenceV.^2);
            
            if nargout > 1 % gradient requested
                grad = zeros(size(x));
                grad(violIdx) = (2 * wV) .* differenceV;
                
                if nargout > 2 % Hessian mult. function requested
                    hessMultVecorFunc = @hessMult;
                end
            end
            
            function hessV = hessMult(v)
                hessV = zeros(size(v));
                hessV(violIdx) = (2 * wV);
            end
                
                    
        end
            
    end
end
