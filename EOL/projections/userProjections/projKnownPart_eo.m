classdef projKnownPart_eo < Projection_eo
    
    
    % Copyright 2011 Eli Osherovich.
    
    
    properties (Access = private)
        KnownMask = []   % A logical mask indicating location of the known part
        KnownPart = []   % Known part values
        LowerMask = []   % A logical mask indicating location of the known lower bound
        LowerIdx = []    % Indices that corresponds to the 'true' values of LowerMask
        LowerBound = []  % Lower bound values
        UpperMask = []   % A logical mask indicating location of the known upper bound
        UpperIdx = []    % Indices that corresponds to the 'true' values of UpperMask
        UpperBound = []  % Upper bound values
    end
    
    
    methods
        
        function self = projKnownPart_eo(knownMask, knownPart, ...
                lowerMask, lowerBound, upperMask, upperBound)
         
            
            % Sanity checks.
            error(nargchk(2, 6, nargin));
            
            validateattributes(knownMask, {'numeric', 'logical'}, {'binary'});
            validateattributes(knownPart, {'numeric'}, {'nonnan', 'finite'});
            assert(nnz(knownMask) == numel(knownPart));
            self.KnownMask = knownMask(:);
            self.KnownPart = knownPart(:);
            
            if nargin > 2
                validateattributes(lowerMask, {'numeric','logical'}, {'binary'});
                validateattributes(lowerBound, {'numeric'}, {'nonnan', 'finite'});
                assert(nnz(lowerMask) == numel(lowerBound));
                self.LowerMask = lowerMask(:);
                self.LowerBound = lowerBound(:);
                self.LowerIdx = find(lowerMask);
            end
            if nargin > 4
                validateattributes(upperMask, {'numeric','logical'}, {'binary'});
                validateattributes(upperBound, {'numeric'}, {'nonnan', 'finite'});
                assert(nnz(upperMask) == numel(upperBound));
                self.UpperMask = upperMask(:);
                self.UpperBound = upperBound(:);
                self.UpperIdx = find(upperMask);
            end
            
            % Do not allow any intersection with the KnownMask.
            if ~isempty(self.KnownMask)
                if ~isempty(self.LowerMask) && any(self.LowerMask & self.KnownMask)
                    error('KnownMask and LowerMask must be disjoint');
                end
                if ~isempty(self.UpperMask) && any(self.UpperMask & self.KnownMask)
                    error('KnownMask and UpperMask must be disjoint');
                end
            end
            
            % If both the LowerMask and UpperMask are provided make sure
            % they pose consistent requirements.
            if ~isempty(self.UpperMask) && ~isempty(self.LowerMask)
                tmpLower = -Inf(size(self.LowerMask));
                tmpLower(self.LowerMask) = self.LowerBound;
                tmpUpper = Inf(size(self.UpperMask));
                tmpUpper(self.UpperMask) = self.UpperBound;
                assert(all(tmpUpper(:) >= tmpLower(:)));
            end
        end
        
        function [xNew, dist, xProj] = doProjection(self, x)
            
            % Sanity check
            assert(numel(x) == numel(self.KnownMask));
            
            % Generate the projection:
            
            % Enforce the known part.
            xNew = x;
            xNew(self.KnownMask) = self.KnownPart;
            
            % Enforce the lower bound.
            if ~isempty(self.LowerMask)
                assert(numel(x) == numel(self.LowerMask));
                violIdx = xNew(self.LowerMask) < self.LowerBound;
                xNew(self.LowerIdx(violIdx)) = self.LowerBound(violIdx);
            end
            
            % Enforce the upper bound.
            if ~isempty(self.UpperMask)
                assert(numel(x) == numel(self.UpperMask));
                violIdx = xNew(self.UpperMask) > self.UpperBound;
                xNew(self.UpperIdx(violIdx)) = self.UpperBound(violIdx);
            end
            
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
