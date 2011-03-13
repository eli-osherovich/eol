function db = calcSNR_db(x_true, x_noisy)
    % CALCSNR_DB - calculate Signal-To-Noise ratio in decibels.
    
    
    % Copyright 2009-2011 Eli Osherovich.
    
    
    noise = x_noisy(:) - x_true(:);
    db = 20*log10(norm(x_true(:))/norm(noise));
    
