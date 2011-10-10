function F = genFFTMatrix (xSize)
    % GENFFTMATRIX - generate FFT matrix F such that F*x(:) = fftn(x)(:).
    
    
    % Copyright 2008-2011 Eli Osherovich.
    
    
    % Some sanity checks.
    validateattributes(xSize, {'numeric'}, ...
        {'integer', 'positive', 'vector'});
    
    F = 1;
    for i = 1:length(xSize)
        F = kron(fft(eye(xSize(i))), F);
    end
