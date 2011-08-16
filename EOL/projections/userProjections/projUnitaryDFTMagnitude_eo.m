classdef projUnitaryDFTMagnitude_eo < Projection_eo
    
    
    
    % Copyright 2011 Eli Osherovich.
    
    
    
    properties (Access = private)
        FMag              % Fourier trasform magnitude
        Fop               % Unitary DFT operator
        ComplexVarsFlag   % An indicator whether x is complex or real
    end
    
    methods
        
        function self = projUnitaryDFTMagnitude_eo (shape, fmag, complexVarsFlag)
            
            % Perform some sanity checks:
            
            % Shape must be a vector of positive integers.
            validateattributes(shape, {'numeric'}, ...
                {'integer', 'positive', 'vector'})
            
            % Fourier magnitude fmag must be an array of nonnegative real
            % values. The number of elements in fmag must match the shape.
            validateattributes(fmag, ...
                {'numeric'}, {'real', 'nonnegative', 'nonnan', 'finite'});
            
            % complexVarsFlag must be a logical scalar.
            validateattributes(complexVarsFlag, ...
                {'logical'}, {'scalar'});
            
            
            self.FMag = fmag(:);
            self.Fop = UnitaryDFT_eo(shape);
            self.ComplexVarsFlag = complexVarsFlag;
        end
        
        function [xNew, dist, xProj] = doProjection(self, x)
            % Introduce local variables instead of object properties.
            % Current MATLAB version does not use JIT for calculations
            % involving object properties.
            fmag = self.FMag;
            F = self.Fop;
            complexVarsFlag = self.ComplexVarsFlag;
            
            % Trnasform x into the Fourier domain.
            xF = F*x;
            
            % Set the correct Fourier magnitue.
            xProj = fmag .* exp(1i*angle(xF));
            
            % Calculate the distance in the projection (Fourier) domain (if
            % requrired).
            if nargout > 1
                dist = norm(xF - xProj);
            end
            
            
            % Convert back to the object domain. We use the fact that F is
            % orthogonal, i.e., inv(F) = F'.
            if complexVarsFlag
                xNew = F'*xProj;
            else
                xNew = real(F'*xProj);
            end
        end
        
        
    end
end
