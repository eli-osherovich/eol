classdef pfComplexAbsEQ_eo < PenaltyFunc_eo
    % COMPLEX_ABS_EQ - Squared L2 norm of (|z| - R).
    % COMPLEX_ABS_EQ computes wieghted (optional) squared L2 norm of (|z| - R).
    
    
    
    % Copyright 2008-2010 Eli Osherovich.
    
    
    
    properties (Access=private)
        r;      % given absolute value
        w = 1;  % weights
    end
    
    methods
        function self = pfComplexAbsEQ_eo(r, w)
            
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
            
            % Calculate z's modulus and phase
            % zModulus = abs(z);
            % zPhase = angle(z);
            [zPhase, zModulus] = cmplx2polC_eo(z);
            
            % Calculate function value.
            val = sum(W .* (zModulus - R).^2);
            
            if nargout > 1, % gradient requested
                %zNormalized = complex(cos(zPhase), sin(zPhase));
                zNormalized = pol2unitcmplxC_eo(zPhase);
                grad = 2*W .* (z - R .* zNormalized);
                
                if nargout > 2 % Hessian mult. function is requested
                    hessMultVecorFunc = @hessMult;
                end
            end
            
            function hessV = hessMult(v)
                rzModRatio = R./(zModulus+eps);
                hessV = W .* (...
                    (2 - rzModRatio) .* v + ...
                    (rzModRatio .* zNormalized.^2) .* conj(v));
            end
        end
    end
end















