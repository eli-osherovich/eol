function [coordsObj, coordsFourier] = gen_DFT_coords(nsamples, phys_size)
%GEN_DFT_COORDS - generates physical coordinates.
% [COORDSOBJ, COORDSFOURIER] = GEN_DFT_COORDS(NSAMPLES, PHYS_SIZE) -
% generates object and Fourier space coordinates that correspond to the
% physical size PHYS_SIZE and to the number of samples NSAMPLES.
%
% The physical object is assumed to be N-dimensional.
%
% The NSAMPLES vector is of length N and i-th entry defines the number of
% samples in the i-th dimension.
%
% PHYS_SIZE is either a scalar or vector of length N. In the former case
% the object is assumed to be a N-dimensional cube of size PHYS_SIZE. If
% PHYS_SIZE is a vector its i-th entry sepcifies the extent of the object
% in the i-th dimension.
%
% The coordinates COORDSOBJ and COORDSFOURIER are returned as a cell arrays
% of length N. Note that the origin is put in the middle, i.e., object's
% extent is given by -PHYS_SIZE/2 : PHYS_SIZE/2
%


% (c) 2010 Eli Osherovich.

    % check inputs validity
    assert (isvector(nsamples))
    assert (isvector(phys_size) || isscalar(phys_size))
    
    if ~isscalar(phys_size)
        assert(length(phys_size) == length(nsamples))
    end
    
    % preallocate space
    coordsObj = cell(length(nsamples), 1);
    coordsFourier = coordsObj;
    
    % generate coordinates with zero at the middle
    for d = 1:length(nsamples)
        n = nsamples(d);
        
        % physical step size
        if ~isscalar(phys_size)
            dx = phys_size(d)/n;
            df = 1/phys_size(d);
        else
            dx = phys_size/n;
            df = 1/phys_size;
        end
        
        % generate coordinates
        tmp = floor(-(n-1)/2) : floor((n-1)/2);
        
        coordsObj{d} = tmp*dx;
        coordsFourier{d} = tmp*df;
    end
    
