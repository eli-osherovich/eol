classdef pfComplexAbsLEQ_eo < PenaltyFunc_eo
    % COMPLEX_ABS_LEQ - Squared L2 norm of (|z| - R) where |z| > R.
    % COMPLEX_ABS_EQ computes wieghted (optional) squared L2 norm of (|z| - R) where |z| > R,
    % thus, forcing |z| to be less than r.
    
    
    % Copyright 2008-2011 Eli Osherovich.

    
    
    properties (Access=private)
        r;      % given absolute value
        w = 1;  % weights
    end
    
    methods
        function self = pfComplexAbsLEQ_eo(r, w)
            
            % Check that R is nonnegative.
            validateattributes(r, {'numeric'}, {'real', 'nonnegative'});
            self.r = r(:);
            
            % Set the weight matrix W (if provided)
            if 2 == nargin
                validateattributes(w, {'numeric'}, {'real', 'nonnegative'});
                self.w = w(:);
            end
        end
    
        function [val, grad, hessMultVectorFunc] = doCalculations(self, z)
            % Introduce local variables instead of object's members.
            % Current MATLAB version does not use JIT for calculations
            % involving object's members.
            R = self.r;
            W = self.w;
            
            % Calculate z's modulus and phase.
            % zModulus = abs(z);
            % zPhase = angle(z);
            [zPhase, zModulus] = cmplx2polC_eo(z);
            
            difference = zModulus - R;

            % Find violations.
            violIdx = difference > 0;

            % Violating values.
            zv = z(violIdx);
            zvModulus = zModulus(violIdx);
            zvPhase = zPhase(violIdx);
            Rv = R(violIdx);
            if isscalar(W)
                Wv = W;
            else
                Wv = W(violIdx);
            end

            % Calculate function value.
            val = sum(Wv .* difference(violIdx).^2);

            if nargout > 1, % gradient requested
                %zvNormalized = complex(cos(zvPhase), sin(zvPhase));
                zvNormalized = pol2unitcmplxC_eo(zvPhase);
                grad = zeros(size(z));
                grad(violIdx) = 2*Wv .* (zv -  Rv.*zvNormalized);
                
                if nargout > 2 % Hessian mult. function is requested
                    hessMultVectorFunc = @hessMult;
                end
            end
            
            function hessV = hessMult(v)
                rzvModRatio = Rv./(zvModulus+eps);
                hessV = zeros(numel(v), 1);
                vv = v(violIdx);
                hessV(violIdx) = Wv .* ( ...
                    (2 - rzvModRatio) .* vv + ...
                    (rzvModRatio.*zvNormalized.^2) .* conj(vv));
            end
        end
    end
end
