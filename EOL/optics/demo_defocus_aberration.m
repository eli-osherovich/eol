function demo_defocus_aberration()
% A defocus demo.

% (c) 2010 Eli Osherovich.
    
    %% set dimensions
    
    % wave length (m)
    lambda = 650e-9; % 650 nm
    
    % image size L x L (m)
    L = 1e-3; % 1 mm
    
    % lens raidus (m)
    R = 1e-2; % 1 cm
    
    % distance to image (when focused) (m)
    zi = 0.1; % 10 cm
    
    % distance to image (actual - defocused) (m)
    za = zi - zi*0.01; % 1 percent defocus 
    
    
    % load image
    I = sqrt(im2double(imread('lena_512.png')));
    
    % generate coordinates 
    [~, c_fourier] = gen_DFT_coords(size(I), L);
    [Fy, Fx] = ndgrid(c_fourier{1}, c_fourier{2});
    
    % lens aperture 
    H = zeros(size(I));
    H((Fx.^2 + Fy.^2) * (lambda*zi)^2 <= R^2) = 1;
    
    % add defocus aberrations to the pupil plane
    W = exp(1i * lambda*zi * (1 - zi/za) *(Fx.^2 + Fy.^2));
    H = H.*W;
    
    
    % generate final image
    If = ifft2(fft2(I).*ifftshift(H));
    
    figure;
    imshow(abs(H));
    
    figure;
    imshow(I.^2);
    
    figure;
    imshow(abs(If).^2);
    
   
    
    
    
