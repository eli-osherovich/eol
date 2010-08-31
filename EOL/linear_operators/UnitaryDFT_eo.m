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
% (vectors) are stored in the properties RANGESIZE and IMAGESIZE. The
% mode of operations is governed by the state of ADJOINTFLAG: false =>
% forward operator, true => adjoint operator. To allow MATLAB's notation
% A*x the MTIMES subroutine is overloaded. To allow the convenient
% notation A', CTRANSPOSE subroutine is overloaded.
    
    
    
% Copyright 2008-2010 Eli Osherovich.
    
    
    
    
    properties (Access = private)
        NormConstFwd;
        NormConstAdj;
        ValidIdx = cell(0);
    end
    
    methods
        function self = UnitaryDFT_eo(rangeSize, imageSize)
            if nargin == 2
                linopArgs = {rangeSize, imageSize};
            elseif nargin == 1
                linopArgs = {rangeSize};
            else
                error('EOL:UnitaryDFT:WrongArgNum', ...
                    'You must provide either one or two arguments');
            end
            self = self@LinearOp_eo(linopArgs{:});
            % We currently allow only zero padding (image size is greater
            % than or equal to the range size). No "chopping" is allowed.
            if any(self.ImageSize < self.RangeSize)
                error('EOL:UnitaryDFT:WrongSize', ...
                      'Range size is greater than the image size');
            end
            
            % Set normalization constant. By multiplying by this constant
            % we make the operator unitary.
            self.NormConstAdj = sqrt(prod(self.ImageSize));
            self.NormConstFwd = 1/self.NormConstAdj;
            
            
            % Set valid index.
            if any(self.ImageSize ~= self.RangeSize)
                self.ValidIdx = cell(1, length(self.RangeSize));
                for k = 1:length(self.ValidIdx)
                    self.ValidIdx{k} = 1:self.RangeSize(k);
                end
            end
        end
        
        function Ax = ApplyAdjoint(self, x)
            x_tmp  = reshape(x, self.ImageSize);
            Ax = ifftn(x_tmp);
            if ~isempty(self.ValidIdx)
                Ax = Ax(self.ValidIdx{:});
            end
            Ax = Ax * self.NormConstAdj;
            Ax = Ax(:);
        end
        
        function Ax = ApplyForward(self, x)
            x_tmp = reshape(x, self.RangeSize);
            Ax = fftn(x_tmp, self.ImageSize);
            Ax = Ax*self.NormConstFwd;
            Ax = Ax(:);
        end
    end
end
