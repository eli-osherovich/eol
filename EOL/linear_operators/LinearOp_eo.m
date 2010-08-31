classdef LinearOp_eo
    
    properties (Access = protected)
        RangeSize
        RangeNumel
        
        ImageSize
        ImageNumel
        
        AdjointFlag = false
    end
    
    methods        
        function self = LinearOp_eo(rangeSize, imageSize)
            if nargin < 1
                error('EOL:LinearOp:WrongInput', ['Not enough input ' ...
                    'arguments']);
            end
            
            % Inputs must be integer positive vectors.
            validateattributes(rangeSize, {'numeric'}, {'integer', ...
                'positive', 'vector'});
            self.RangeSize = rangeSize;
            
            % Set image shape (size).
            if nargin > 1
                self.ImageSize = imageSize;
            else
                self.ImageSize = rangeSize;
            end
            
            self.RangeNumel = prod(self.RangeSize);
            self.ImageNumel = prod(self.ImageSize);
        end
        
        % Toggle the conjugate transpose flag on a call to ctranspose. By
        % defining this function we allow correct treatment of the '
        % (apostrophe) operator. The operator works correctly on an array
        % of LinearOp_eo(s).
        function self = ctranspose(self)
            [self.AdjointFlag] = deal(~[self.AdjointFlag]);
            self = self.';
        end
        
        function Ax = mtimes(self, x)
            if isscalar(self)
                if self.AdjointFlag
                    Ax = ApplyAdjoint(self, x);
                else
                    Ax = ApplyForward(self, x);
                end
            elseif ndims(self) == 2
                Ax = [];
                for i = 1:size(self, 1)
                    Ax_row = 0;
                    for j = 1:size(self, 2)
                        if self(i,j).AdjointFlag
                            Ax_row = Ax_row +  ApplyAdjoint(self(i,j), x);
                        else
                            Ax_row = Ax_row + ApplyForward(self(i,j), x);
                        end
                    end
                    Ax = [Ax; Ax_row];
                end
            else
                error('Oooops');
            end
        end
   
        
        % To concatenate linear operators horizonally we must check that the
        % image dimesnions of all arguments is the same. There are two
        % possibilities: first, to check the shape (size) of the image;
        % second, to check only that the number of elements in the image is
        % the same. Here we use the second approach.
        function self = horzcat(varargin)
            self = builtin('horzcat', varargin{:});
            % Make sure that the image is compatible.
            if ~all(self(1).ImageNumel == [self.ImageNumel])
                error('EOL:LinearOp:horzcat:dimensionMismatch', ...
                    'CAT arguments dimensions are not consistent');
            end
        end
        
        % Similarly to the horizontal concatenation, the vertical
        % concatenation must check that the range is compatible.
        function self = vertcat(varargin)
            self = builtin('vertcat', varargin{:});
            % Make sure that the range is compatible.
            if ~all(self(1).RangeNumel == [self.RangeNumel])
                error('EOL:LinearOp:vertcat:dimensionMismatch', ...
                    'CAT arguments dimensions are not consistent');
            end
        end
        
        % General cat function: do not allow creation of ND arrays. 
        % It must be a matrix.
        function self = cat(varargin)
            switch varargin{1}
                case 1
                    self = vertcat(varargin{2:end});
                case 2
                    self = horzcat(varargin{2:end});
                otherwise
                    error('EOL:LinearOp:cat:dimensionError', ...
                        'CAT can create matrices only. Higher dimensional arraya are not allowed');
            end
        end
    end
    
    methods (Abstract)
        Ax = ApplyAdjoint(self, x)
        Ax = ApplyForward(self, x)
    end
end
