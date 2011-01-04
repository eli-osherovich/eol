classdef KnownPart_eo < AffineOp_eo
    
    % Usage example:
    % --------------
    % a = magic(5)
    % m = false(5);
    % m(1:2, 1:2) = true;
    % KP = KnownPart_eo(m, a(m));
    % x = rand(25-4, 1);
    % reshape(KP*x, size(a))
    % KP'*a
    
    % Inputs:
    % -------
    % To create an object you have to specify a mask (logical array)
    % whose non-zero (true) values corresponding the locations of the known
    % part. In addition the known values must be provided.
    %
    % Output:
    % -------
    % Output is always a column vector.
    %
    
    
    
    % Copyright 2010 Eli Osherovich.
    
    
    properties
        KnownMask
        KnownValues
    end
    
    methods
        function self = KnownPart_eo(knownMask, knownValues)
            
            % Zeros in knownMask are locations of the unknown signal part.
            % While non-zeros entries of the mask are locations of the
            % known part. The corresponding values are given in
            % knownValues.
            
            % First, create the Support operator for known values.
            S = Support_eo(~knownMask);
            
            % Make sure that the additive part b is of appropriate size.
            b = zeros(size(knownMask));
            b(knownMask) = knownValues;
            self = self@AffineOp_eo(S', b);
            
        end
    end
end
