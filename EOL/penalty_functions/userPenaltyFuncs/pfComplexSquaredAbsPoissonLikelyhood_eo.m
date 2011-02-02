classdef pfComplexSquaredAbsPoissonLikelyhood_eo < PenaltyFunc_eo
   
    
    
    % Copyright 2008-2011 Eli Osherovich.
    
    
    
    properties (Access=private)
        r;      % given absolute value
        w = 1;  % weights
    end
    
    methods
        function self = pfComplexSquaredAbsPoissonLikelyhood_eo(r, w)
            
            % Check that R is nonnegative.
            validateattributes(r, {'numeric'}, {'real', 'nonnegative'});
            self.r = r(:);
            
            % Set the weight matrix W (if provided)
            if 2 == nargin
                validateattributes(w, {'numeric'}, {'real', 'nonnegative'});
                self.w = w(:);
            end
        end
    
        function [val, grad, hessMultVecorFunc] = doCalculations(self, z)
            % Introduce local variables instead of object's members.
            % Current MATLAB version does not use JIT for calculations
            % involving object's members.
            R = self.r;
            W = self.w;
            
            % calculate |z|^2
            zSquaredModulus = z.*conj(z);
                        
            % Calculate function value.
            val = sum(W .* (zSquaredModulus(:)-R(:).*log(zSquaredModulus(:))));
            
            if nargout > 1, % gradient requested
                
                grad = 2*W .* (z - R./conj(z));
                
                if nargout > 2 % Hessian mult. function is requested
                    hessMultVecorFunc = @hessMult;
                end
            end
            
            function hessV = hessMult(v)
                hessV = 2*W .* (v + R./conj(z.^2).*conj(v));
            end
        end
    end
end















