classdef CircularConv_eo < LinearOp_eo
    % A circular N-dimensional convolution class.
    
    
    
    % Copyright 2009-2011 Eli Osherovich.
    
    
    
    properties
        Kernl
        KernlFourier
        ImageShape
    end
    
    methods 
        function self = CircularConv_eo(kernl, imageShape)
            self = self@LinearOp_eo(prod(imageShape), prod(imageShape));
            
            % Save the original kernel (why?) and its Fourier trasform (of
            % appropriate shape).
            %
            % TODO: check whether the kernel's separability can sepeed up
            % the calculations.
            self.Kernl = kernl;
            
            % For correct application, the kernel and image must be of same
            % rank. Hence, we add singletone dimenstions if their ranks are
            % different.
            kernlRank = ndims(kernl);
            imageRank = length(imageShape);
            
            kernlShape = [size(kernl), ones(1, imageRank - kernlRank)];
            imageShape = [imageShape, ones(1, kernlRank - imageRank)];
            
            self.ImageShape = imageShape;
            
            % Verify that the kernel size is not greater than image size.
            if any(kernlShape > imageShape)
                error('EOL:CircularConv:CircularConv:WrongArg', ...
                    'Kernel size must not exceed image size');
            end
            
            % Center of the kernel.
            kernlCenter = floor(kernlShape/2);

            % Pad kernel with zerros (if necessary) and shift so that the
            % center will be at the origin. 
            kernlPadded = padarray(reshape(kernl, kernlShape), ...
                imageShape - kernlShape, 0, 'post');
            kernlPadded = circshift(kernlPadded, -kernlCenter);
            
            % Save kernel's Fourier transform.
            self.KernlFourier = fftn(kernlPadded);
        end
        
        function Ax = ApplyForward(self, x)
            x_tmp = reshape(x, self.ImageShape);
            Ax = ifftn(self.KernlFourier .* fftn(x_tmp));
        end
        
        function Ax = ApplyAdjoint(self, x)
            x_tmp = reshape(x, self.ImageShape);
            Ax = ifftn(conj(self.KernlFourier) .* fftn(x_tmp));
        end
    end
end
