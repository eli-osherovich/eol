classdef pfComplexSquaredAbsEQ_eo < PenaltyFunc_eo
    % COMPLEX_ABS_SQUARED_EQ - penalty function that forces squared z modulus
    % to be equal to R, i.e., |z|^2 = R
    
    
    
    % Copyright 2008-2011 Eli Osherovich.
    
    properties
        r;           % given absolute value squared
        w = 1;       % weights.
    end
    
    methods
        function self = pfComplexSquaredAbsEQ_eo(r, w)
            
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
            
            % Calculate z's abs value squared.
            zAbsSquared = z.*conj(z);
            assert(isreal(zAbsSquared));
            
            % Calculate function value.
            val = sum(W .* (zAbsSquared-R).^2);
            
            
            if nargout > 1, % gradient requested
                grad = 4*W .* z.*(zAbsSquared - R);
                
                if nargout > 2 % Hessian mult. function is requested.
                    hessMultVectorFunc = @hessMult;
                end
            end
            
            function hessV = hessMult(v)
                hessV = 4*W .* (...
                    (2*zSquaredModulus - R) .* v + ...
                    z.^2 .* conj(V));
            end
        end
    end
end
