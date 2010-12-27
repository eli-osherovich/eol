classdef KnownPart_eo < LinearOp_eo
    % Affine operator that implements Ax+b where A is a padding operator
    % that creates some "free" space for the known part (given by b);
    %
    % Usage example:
    % --------------
    % a = magic(5)
    % m = false(5);
    % m(1:2, 1:2) = true;
    % KV = KnownPart_eo(m, a(m));
    % x = rand(25-4, 1);
    % reshape(KV*x, size(a))
    % KV'*a
    
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
    
    
    properties (Access = private)
        KnownMask
        KnownValues
    end
    
    methods
        function self = KnownPart_eo(knownMask, knownValues)
            validateattributes(knownMask, {'logical'}, {'nonempty'});
            
            % Zeros in knownMask are locations of the unknown signal part.
            % While non-zeros entries of the mask are locations of the
            % known part. The corresponding values are given in
            % knownValues.
            self = self@LinearOp_eo(nnz(~knownMask), numel(knownMask));
            
            % Save the known values and their locations.
            self.KnownMask = knownMask(:);
            self.KnownValues = knownValues(:);
        end
        
        % Forward opterator: pad x with zeros and fill those new locations
        % with known values.
        function Ax = ApplyForward(self, x)
            Ax = zeros(numel(self.KnownMask), 1);
            Ax(self.KnownMask) = self.KnownValues;
            Ax(~self.KnownMask) = x;
        end
        
        % Adjoint operator: extract only unknown entries (variables).
        function Ax = ApplyAdjoint(self, x)
            Ax = x(~self.KnownMask);
        end
    end
end
