classdef pfProjGEQL2Squared_eo < PenaltyFunc_eo
    
    
    
    % Copyright 2008-2010  Eli Osherovich.
    
    
    properties (Access=private)
        dir = 1
        threshold = 0
        w = 1
    end
    
    methods
        function self = pfProjGEQL2Squared_eo(dir, thr, w)
            switch nargin
                case 0
                    % Do nothing
                    % This is the default case: the real part of x
                    % projected onto the non-negative orthant.
                case 1
                    validateattributes(dir, {'numeric'}, {'finite'})
                    self.dir = dir(:);
                case 2
                    validateattributes(dir, {'numeric'}, {'finite'});
                    validateattributes(thr, {'numeric'}, {'real'});
                    self.dir = dir(:);
                    self.threshold = thr(:);
                case 3
                    validateattributes(dir, {'numeric'}, {'finite'});
                    validateattributes(thr, {'numeric'}, {'real'});
                    validateattributes(w, {'numeric'}, {'real', 'nonnegative'});
                    self.dir = dir(:);
                    self.threshold = thr(:);
                    self.w = w(:);
            end
            % Make sure that the norm of each entry in DIR is unity.
            self.dir = self.dir ./ abs(self.dir);
            
        end

        function [val, grad, hessMultVecorFunc] = doCalculations(self, x)
            Dir = self.dir;
            Thr = self.threshold;
            W = self.w;
           
            % Violations: Projection of X on DIR < threshold
            difference = real(Dir).*real(x) + imag(Dir).*imag(x) - Thr;
            violIdx = difference < 0;
            
            % Parameters' values at violating indices. 
            if isscalar(Dir)
                dirV = Dir;
            else
                dirV = Dir(violIdx);
            end
            if isscalar(W)
                wV = W;
            else
                wV = W(violIdx);
            end
            differenceV = difference(violIdx);
            wVDiff = wV .* differenceV;
            
            val = sum( wVDiff .* differenceV);
            
            if nargout > 1 % gradient requested
                grad = zeros(size(x));
                grad(violIdx) = 2 * wVDiff .* dirV;
                if nargout > 2 % Hessian mult. function requested
                    hessMultVecorFunc = @hessMult;
                end
            end
            
            function hessV = hessMult(v)
                hessV = zeros(size(v));
                hessV(violIdx) = wV .* (v(violIdx) + dirV.^2 .* conj(v(violIdx)));
            end
                
                    
        end
            
    end
end
