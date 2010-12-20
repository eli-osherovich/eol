classdef Support_eo < LinearOp_eo
    % Linear operator that implements support (location that are non-zeros).
    %
    % Usage example:
    % --------------
    % a = magic(5)
    % m = false(5);
    % m(1:2, 1:2) = true;
    % S = Support_eo(m);
    %
    % Truncation:
    % t = S*a;
    % 
    % Zero padding:
    % p = S'*t;
    
    % Inputs:
    % -------
    % To create a support object you have to specify a mask (logical array)
    % whose non-zero (true) values will indicate the support (known
    % locations of a signal's non-zero values).
    %
    % Output:
    % -------
    % Output is always a column vector.
    %
    
    
    
    % Copyright 2010 Eli Osherovich.
    
    
    properties (Access = private)
        SupportMask
    end
    
    methods
        function self = Support_eo(mask)
            validateattributes(mask, {'logical'}, {'nonempty'});
            
            % The operator is applied to signals whose shape is equal to
            % that of the mask. Hence, the range dimensionality is given by
            % the total number of elements in the mask. While image's
            % dimensionality is the number of non-zeros in the mask.
            self = self@LinearOp_eo(numel(mask), nnz(mask));
            
            % Save the support mask.
            self.SupportMask = mask(:);
        end
        
        % Forward opterator: extract x's entries defined by the mask's
        % non-zeros.
        function Ax = ApplyForward(self, x)
            Ax = x(self.SupportMask);
        end
        
        % Adjoint operator: distribute x's entries into the locations
        % defined by non-zeros of the mask.
        function Ax = ApplyAdjoint(self, x)
            Ax = zeros(self.RangeNumelOriginal, 1);
            Ax(self.SupportMask) = x;
        end
    end
end
