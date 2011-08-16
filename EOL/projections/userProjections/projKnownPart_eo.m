classdef projKnownPart_eo < Projection_eo
    
    
    % Copyright 2011 Eli Osherovich.
    
    
    properties (Access = private)
        KnownMask   % A logical mask indicating location of the known part
        KnownPart   % Known part values
    end
    
    
    methods
        
        function self = projKnownPart_eo(knownMask, knownPart)
         
            % Sanity checks.
            validateattributes(knownMask, {'logical'}, {'nonempty'});
            validateattributes(knownPart, {'numeric'}, {'nonnan', 'finite'});
            
            assert(nnz(knownMask) == numel(knownPart));
            
            self.KnownMask = knownMask(:);
            self.KnownPart = knownPart(:);
        end
        
        function [xNew, dist, xProj] = doProjection(self, x)
            
            % Sanity check
            assert(numel(x) == numel(self.KnownMask));
            
            % Generate the projection.
            xNew = x;
            xNew(self.KnownMask) = self.KnownPart;
            
            % Calculate the rest of output argumetns.
            if nargout > 1
                dist = norm(xNew(:) - x(:));
                
                if nargout > 2
                    xProj = xNew;
                end
            end
        end
    end
end
