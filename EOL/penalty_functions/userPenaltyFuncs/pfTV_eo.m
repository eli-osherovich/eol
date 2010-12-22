classdef pfTV_eo < PenaltyFunc_eo
    % TV - total variation
    
    
    % Copyright 2008-2010 Eli Osherovich.
    
    
    properties 
        shape % shape of the signal
        nDims % number of dimensions (rank) of the signal
        epsilon = 0.001 % smoothing parameter used in |x| approximation.
    end
    
    methods
        function self = pfTV_eo(shape, e)
            % Make sure that the shape is an integer vector.
            validateattributes(shape, {'numeric'}, {'integer', 'vector'});
            self.shape = shape;
            self.nDims = length(shape);
            
            % Set smoothing parameter (epsilon) if provided.
            if 2 == nargin
                validateattributes(e, {'numeric'}, {'real', 'positive'});
                self.epsilon = e;
            end
        end
        
        function [val, grad, hessMultVecorFunc] = doCalculations(self, x)
            % Reshape x to its original (correct) shape.
            x = reshape(x, self.shape);
            
            % Preallocate space for derivatives of x
            % (i'th column contains the derivative in the i'th direction).
            if isreal(x)
                Dx = zeros(numel(x), self.nDims);
            else
                Dx = complex(zeros(numel(x), self.nDims), zeros(numel(x), nDims));
            end
            
            % Calculate derivatives in all directions.
            for d = 1:self.nDims
                tmp = Dop(x, d);
                Dx(:, d) = tmp(:);
            end
            
            % Smoothed approximation to abs(grad(x))
            abs_gradX = sqrt(sum(abs(Dx).^2,2) + self.epsilon);
            
            % A simple "normalization": make sure the value vanishes for a
            % constant input.
            val = sum(abs_gradX - sqrt(self.epsilon));
            
            if nargout > 1, % gradient requested
                grad = zeros(self.shape);
                
                for d = 1:self.nDims
                    tmp = Dx(:, d)./abs_gradX;
                    tmp = reshape(tmp, self.shape);
                    grad = grad + DopAdj(tmp, d);
                end
                
                % Finally reshape grad to a column vector.
                grad = grad(:);
                
                if nargout > 2 % Hessian mult. function is requested
                    hessMultVecorFunc = @hessMult;
                end
            end
        end
        
        function hessV = hessMult(v)
            error('Not implemented : need a review for complex case.');
            
            if isnumeric(V) % one vector
                assert(isvector(V));
                
                % pre-allocate space
                hessV = zeros(size(V));
                
            elseif iscell(V) % multiple vectors
                % pre-allocate space
                hessV = cell(size(V));
                
                for s = 1:size(V, 2)
                    for i = 1:realDims
                        tmp = zeros(size(x));
                        for j = 1:realDims
                            Ds = Dop(reshape(V(:,s), size(x)), j);
                            if i==j
                                tmp(:) = tmp(:) + (abs_gradX.^2 - Dx(:,i).^2)./abs_gradX.^3.*Ds(:);
                            else
                                tmp(:) = tmp(:) - Dx(:,i).*Dx(:,j)./abs_gradX.^3.*Ds(:);
                            end
                        end
                        tmp = DopAdj(tmp, i);
                        hessV(:,s) = hessV(:,s) + tmp(:);
                    end
                end
            else
                error('Wrong dimensionality')
            end
        end
    end
end



function y = Dop(x, curr_dim)
% DOP - Differention operator
% Y = DOP(X, CUR_DIM) - performs differentiation of X along dimension
% CUR_DIM.
% The derivative is approximated by forward difference. 
% At the end boundary the derivative is assumed to be zero (corresponds to
% mirror boundary conditions).
%
% Example: given a vector x = [x1, x2, x3, x4]' after application of the
% operator one gets [x2-x1, x3-x2, x4-x3, 0]';

x_size = size(x);

% pre allocate space according the type of x (real or complex)
if isreal(x)
    y = zeros(size(x));
else
    y = complex(zeros(size(x)), zeros(size(x)));
end

% calculate valid indices
% all indices are valid except the last index in the direction of
% differentiation
idx = cell(1, length(x_size));
for d = 1:length(x_size)
    idx{d} = 1:x_size(d);
end
idx{curr_dim} = 1:x_size(curr_dim)-1;

% compute the derivative (forward difference)
y(idx{:}) = diff(x,1, curr_dim);
end


function y = DopAdj (x, curr_dim)
% DOPADJ - Adjoint differentiation operator
% Y = DOPADJ (X, CURR_DIM) - performs adjoint differentiation of X along
% dimension CURR_DIM.
%
% Example: given a vector y=[y1, y2, y3, y4]' after application of the
% operator one gets [-y1, y1-y2, y2-y3, y3];


x_size = size(x);
idx1 = cell(1, length(x_size));
for k = 1:length(x_size)
    idx1{k} = 1:x_size(k);
end
idx_end = idx1;

idx1{curr_dim} = 1;
idx_end{curr_dim} = x_size(curr_dim);

y = cat(curr_dim, -x(idx1{:}), -diff(x, 1, curr_dim));
y(idx_end{:}) = y(idx_end{:}) + x(idx_end{:});
end
