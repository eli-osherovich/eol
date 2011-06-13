classdef pfComplexAbsEQ_eo < PenaltyFunc_eo
    % COMPLEX_ABS_EQ - Squared L2 norm of (|z| - R).
    % COMPLEX_ABS_EQ computes wieghted (optional) squared L2 norm of (|z| - R).
    
    
    
    % Copyright 2008-2011 Eli Osherovich.
    
    
    
    properties (Access=private)
        r;              % given absolute value
        w = 1;          % weights
        phaseToUse = 0; % phase to use when z's magnitude is zero.
    end
    
    methods
        function self = pfComplexAbsEQ_eo(r, ph, w)
            
            % Check that R is nonnegative.
            validateattributes(r, {'numeric'}, {'real', 'nonnegative'});
            self.r = r(:);
            
            % Set phaseToUse (if provided).
            if nargin > 1 && ~isempty(ph)
                if isscalar(ph)
                    validateattributes(ph, {'numeric'}, {'real'});
                else
                    validateattributes(ph(:), {'numeric'}, {'size', [numel(r),1]});
                end
                self.phaseToUse = ph(:);
            end
            
            % Set the weight matrix W (if provided).
            if nargin == 3
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
            % In general z shall not be real. However, Matlab's
            % implementation of the Fourier transform CAN produce real
            % output even when the input was real.
            if isreal(z)
                zModulus = abs(z);
                zPhase = zeros(size(z));
            else
                % zModulus = abs(z);
                % zPhase = angle(z);
                [zPhase, zModulus] = cmplx2polC_eo(z);
            end
           
            
            % The phase is undetermined where zModules = 0. Hence, we have
            % to use some reasonable value. The best approach would
            % probably be to average the phases of the neighbors. Here we
            % use somewhat simpler remedy: use some predefined/given
            % values. Note that by default the phase is set to zero if 
            % |z| = 0.
            PhaseToUse = self.phaseToUse;
            zeroIdx = zModulus == 0;
            if isscalar(PhaseToUse) 
                if PhaseToUse ~= 0
                    zPhase(zeroIdx) = PhaseToUse;
                end
            else
                zPhase(zeroIdx) = PhaseToUse(zeroIdx);
            end
            
            % Calculate function value.
            val = sum(W .* (zModulus - R).^2);
            
            if nargout > 1, % gradient requested
                %zNormalized = complex(cos(zPhase), sin(zPhase));
                zNormalized = pol2unitcmplxC_eo(zPhase);
                grad = 2*W .* (z - R .* zNormalized);
                
                if nargout > 2 % Hessian mult. function is requested
                    hessMultVectorFunc = @hessMult;
                end
            end
            
            function hessV = hessMult(v)
                rzModRatio = R./(zModulus+eps);
                hessV = W .* ( ...
                    (2 - rzModRatio) .* v + ...
                    (rzModRatio .* zNormalized.^2) .* conj(v));
            end
        end
    end
end















