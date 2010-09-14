classdef pfComplexAbsEQ_eo < PenaltyFunc_eo
    % COMPLEX_ABS_EQ - Squared L2 norm of (|z| - R).
    % COMPLEX_ABS_EQ computes wieghted (optional) squared L2 norm of (|z| - R).
    
    
    
    % Copyright 2008-2010 Eli Osherovich.
    
    
    
    properties
        r
        w = 1;
    end
    
    methods
        function self = pfComplexAbsEQ_eo(r, w)
            
            % Chech that R is nonnegative.
            validateattributes(r, {'numeric'}, {'real', 'nonnegative'});
            self.r = r(:);
            
            % Set the weight matrix W (if provided)
            if 2 == nargin
                validateattributes(w, {'numeric'}, {'real', 'nonnegative'});
                self.w = w(:);
            end
        end
        
        function [val, grad, hessMultVecorFunc] = doCalculations(self, z)
            % Calculate z's modulus and phase
            % zModulus = abs(z);
            % zPhase = angle(z);
            [zPhase, zModulus] = cmplx2polC_eo(z);
            
            % Calculate function value.
            val = sum(self.w .* (zModulus - self.r).^2);
            
            if nargout > 1, % gradient requested
                %zNormalized = complex(cos(zPhase), sin(zPhase));
                zNormalized = pol2unitcmplxC_eo(zPhase);
                grad = 2*self.w .* (z -  self.r .* zNormalized);
                
                if nargout > 2 % Hessian mult. function is requested
                    hessMultVecorFunc = @hessMult;
                end
            end
            
            function hessV = hessMult(v)
                rzModRatio = self.r./(zModulus+eps);
                hessV = self.w .* (...
                    (2 - rzModRatio) .* v + ...
                    (rzModRatio .* zNormalized.^2) .* conj(v));
            end
        end
    end
end















