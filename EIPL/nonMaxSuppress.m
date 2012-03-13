function maxImSup = nonMaxSuppress(im, neighborMask)
    % NONMAXSUPPRESS - non-maximum suppression for 2D arrays.
    
    
    % Remarsk:
    % --------
    % Im is assumed to be non-negative. 
    
    
    
    
    % Copyright 2012 Eli Osherovich.
    
    
    
    
    % Sanity checks.
    narginchk(2, 2);
    
    validateattributes(im, {'numeric'}, {'nonnegative', '2d'});
    validateattributes(neighborMask, {'numeric'}, {'binary', '2d'});
    
    
    % Maximum is the last element in a sorted array.
    order = nnz(neighborMask);
    maxIm = ordfilt2(im, order, neighborMask);
    maxIdx = maxIm == im;
    
    % Suppress non-maximum values.
    maxImSup = zeros(size(im));
    maxImSup(maxIdx) = im(maxIdx);
    
    
