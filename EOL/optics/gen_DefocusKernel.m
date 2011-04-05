function [h] = gen_DefocusKernel(nsamples, phys_size, zi, za, lambda, R)
%GEN_DEFOCUSKERNEL - generate de-focus kernel.
% H = GEN_DEFOCUSKERNEL(NSAMPLES, PHYS_SIZE, ZI, ZA, LAMBDA) generate
% de-focus kernel H of size  NSAMPLES(1) x NSAMPLES(2). The kernel
% corresponds to an image of physical dimensions given by PHYS_SIZE.
%
% The correct (focused) image should be formed at a distance ZI from exit
% pupil. Instead, it is measured at a distance ZA. 
% 
% LAMBDA denotes the wavelength used. 
%
% R is the aperture radius. This argument is optional. If missing it will
% be assumed that the radius is large enought and does not affect the
% image.
%
% All physical dimensions are given in meters.
%
% Typical values are:
% NSAMPLES = [512, 512];
% PHYS_SIZE = 1e-3; % 1 mm
% ZI = 0.1; % 10 cm
% ZA = ZI - ZI*0.01; % 1 percent focusing error 
% LAMBDA = 650e-9; % 650 nm
% R = 1e-2; % 1 cm

% (c) 2010 Eli Osherovich.



% exit pupil contains (scaled) Fourier transform of the image
% generate physical coordiantes of the exit pupil
[~, c_fourier] = gen_DFT_coords(nsamples, phys_size);

% verify validity of the coordinates
for i = 1:length(c_fourier)
    tmp = ifftshift(c_fourier{i});
    assert(tmp(1) == 0);
end
[Fy, Fx] = ndgrid(c_fourier{1}, c_fourier{2});
 
% aperture image (exit pupil)
H = zeros(nsamples);
H((Fx.^2 + Fy.^2) * (lambda*zi)^2 <= R^2) = 1;
    
% add defocus aberrations to the pupil plane
W = exp(1i * lambda*zi * (1 - zi/za) *(Fx.^2 + Fy.^2));
H = H.*W;


% figure;
% imshow(abs(H));
    
    
% generate defocus kernell
h = fftshift(ifft2(ifftshift(H)));
