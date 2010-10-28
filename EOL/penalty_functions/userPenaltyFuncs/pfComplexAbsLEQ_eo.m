classdef pfComplexAbsLEQ_eo < PenaltyFunc_eo
    % COMPLEX_ABS_LEQ - penalty function that forces z modulus to be less then
    % or equal to R

    
    % Copyright 2008-2010 Eli Osherovich.

    
    properties (Access = private)
        r
        w = 1;
    end
    
    methods
        function self = pfComplexAbsLEQ_eo(r, w)
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
            % Introduce local variables instead of object's members.
            % Current MATLAB version does not use JIT for calculations
            % involving object's members.
            R = self.r;
            W = self.w;
            
            % Calculate z's modulus.
            zModulus = abs(z);
           

            difference = zModulus - R;

            % Find violations.
            violIdx = difference > 0;

            % Violating values.
            zv = z(violIdx);
            zvModulus = zModulus(violIdx);
            rv = R(violIdx);

            val = sum(W.*difference(violIdx).^2);

            if nargout > 1, % gradient requested
                grad = zeros(size(z));
                          
                zvNormalized = exp(1i*angle(zv)); % slow?
                grad(violIdx) = 2*W .* (zv -  rv.*zvNormalized);
                
                if nargout > 2, % Hessian mult. function is requested
                    hessMultVecorFunc = @hessMult;
                end
            end
            
            function hessV = hessMult(v)
                rzvMod_ratio = rv./(zvModulus+eps);
                hessV = zeros(numel(v), 1);
                vv = v(violIdx);
                hessV(violIdx) = W .*( ...
                    (2 - rzvMod_ratio) .*vv + ...
                    (rzvMod_ratio.*zvNormalized.^2) .* conj(vv));
            end
        end
    end
end
