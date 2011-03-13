function snr = calcSNR(x_true, x_noisy)
    % CALCSNR - calculate Signal-To-Noise ratio.
    
    
    % Copyright 2009-2011 Eli Osherovich.
    
    
    
    noise = x_noisy(:) - x_true(:);
    snr = norm(x_true(:))/norm(noise);
    
