classdef UnitaryDFT_eo < LinearOp_eo
% Unitary DFT class.
%
% Usage:
% -------
% without padding:
% x  = rand(5);
% A = UnitaryDFT_eo(size(x));
% Ax = A*x;
% x2 = A'*Ax;
%
% with padding: 
% x = rand(5);
% A =  UnitaryDFT_eo(size(x), 2*size(x));
% Ax = A*x;
% x2 = A'*Ax;   
    
    
% Inputs:
% -------
% To create on object you have to specify either one or two inputs: the
% shape of the original X (range) and the shape of the result AX
% (image). If the second input is omitted then the shape of the image is
% equal to the shape of the range. Note that only zero-padding is
% allowed, i.e., the size of the image must be greater than or equal to
% the size of the range (element-wise).  
%
% To apply the operator, simply use regular MATLAB's notation A*x. To
% create or apply and adjoint operator use A', again, this is the regular
% notation used in MATLAB.
%
%
% Output:
% -------
% Output is always a column vector.
%
%
% Implementation details:
% -----------------------
% The shape (size) of the range and input are mandatory. These values
% (vectors) are stored in the properties RANGESHAPE and IMAGESHAPE. The
% mode of operations is governed by the state of ADJOINTFLAG: false =>
% forward operator, true => adjoint operator. To allow MATLAB's notation
% A*x the MTIMES subroutine is overloaded. To allow the convenient
% notation A', CTRANSPOSE subroutine is overloaded.
    
    
    
% Copyright 2010 Eli Osherovich.
    
    
    
    
    properties (Access = private)
        RangeShape
        ImageShape
        NormConstFwd;
        NormConstAdj;
        ValidIdx = cell(0);
    end
    
    methods
        function self = UnitaryDFT_eo(rangeShape, imageShape)
            if nargin == 2
                validateattributes(rangeShape, {'numeric'}, {'integer', ...
                    'positive', 'vector'});
                validateattributes(imageShape, {'numeric'}, {'integer', ...
                    'positive', 'vector'});
                
                linopArgs = {prod(rangeShape), prod(imageShape)};
                
            elseif nargin == 1
                
                validateattributes(rangeShape, {'numeric'}, {'integer', ...
                    'positive', 'vector'});
               
                linopArgs = {prod(rangeShape)};
                imageShape = rangeShape;
            else
                error('EOL:UnitaryDFT:WrongArgNum', ...
                    'You must provide either one or two arguments');
            end
            self = self@LinearOp_eo(linopArgs{:});
            % We currently allow only zero padding (image size is greater
            % than or equal to the range size). No "chopping" is allowed.
            if any(imageShape < rangeShape)
                error('EOL:UnitaryDFT:WrongSize', ...
                      'Range size is greater than the image size');
            end
            self.RangeShape = rangeShape;
            self.ImageShape = imageShape;
                        
            % Set normalization constant. By multiplying by this constant
            % we make the operator unitary.
            self.NormConstAdj = sqrt(prod(imageShape));
            self.NormConstFwd = 1/self.NormConstAdj;
            
            
            % Set valid index.
            if any(imageShape ~= rangeShape)
                self.ValidIdx = cell(1, length(rangeShape));
                for k = 1:length(self.ValidIdx)
                    self.ValidIdx{k} = 1:rangeShape(k);
                end
            end
        end

        % Forward operator: pad matrix with zeros and apply forward Fourier
        % transform. Note that zero padding is done implicitly by Matlab's
        % fft function (which is, of course, stolen fftw).
        function Ax = ApplyForward(self, x)
            x_tmp = reshape(x, self.RangeShape);
            Ax = fftn(x_tmp, self.ImageShape);
            Ax = Ax*self.NormConstFwd;
            Ax = Ax(:);
        end
        
        % Adjoint operator: apply inverse Fourier transform (which is
        % equivalent to the adjoint) and clip the valid area which is, in
        % general, smaller due to zero padding.
        function Ax = ApplyAdjoint(self, x)
            x_tmp  = reshape(x, self.ImageShape);
            Ax = ifftn(x_tmp);
            if ~isempty(self.ValidIdx)
                Ax = Ax(self.ValidIdx{:});
            end
            Ax = Ax * self.NormConstAdj;
            Ax = Ax(:);
        end
        
        
    end
end
