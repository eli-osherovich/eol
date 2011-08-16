classdef Projection_eo
    
    
    
    % Copyright 2011 Eli Osherovich.
    
    
    % The class is an abstract interface for actual implementations.
    % Each implementatoin must provide 'doProjection' function that accepts
    % current x and returns up to 3 outputs:
    %
    % xNew - the new value of x after the projection (in the x's original
    % domain)
    %
    % dist - distance between the current x and its projection xProj in the
    % domain of projection. 
    %
    % xProj the projection of x in the domain of projection.
    
    
    
    methods (Abstract)
        [xNew, dist, xProj] = doProjection(self, x)
    end
end
