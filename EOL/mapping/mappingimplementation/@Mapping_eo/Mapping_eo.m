classdef Mapping_eo
    
    
    
    % Copyright 2010-2012 Eli Osherovich.
    
    
    
    properties
        RangeNumel
        ImageNumel
        MinusFlag = false
    end
    
    methods
        function self = Mapping_eo (rangeNumel, imageNumel)
            
            % Verify that we have either one or two arguments.
            narginchk(1, 2);
            
            if 1 == nargin
                % Assume imageNumel = rangeNumel if only one argument
                % provided.
                imageNumel = rangeNumel;
            end
            
            % Verify that the arguments are positive integers.
            validateattributes(rangeNumel, {'numeric'}, {'integer', ...
                'positive', 'scalar'});
            validateattributes(imageNumel, {'numeric'}, {'integer', ...
                'positive', 'scalar'});
            
            % Create an object.
            self.RangeNumel = rangeNumel;
            self.ImageNumel = imageNumel;
        end
    end
    
    methods (Abstract)
        Ax = ApplyMapping(self, x);
        Jx = MultJacobian(self, x, xCurrent);
        Jcx = MultConjJacobian(self, x, xCurrent);
    end
end
