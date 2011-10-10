function sizeStr = size2str(sz)
    % SIZE2STR - create a string for given size vector.
    
    % Given a size vector SZ = [N1, N2, N3, ...] it creates the string 
    % 'N1xN2xN3x...'.
    
    
    % Copyright 2011 Eli Osherovich.
    
    
    % Some sanity checks.
    validateattributes(sz, {'numeric'}, {'integer', 'nonnegative', 'vector'});
    
    
    sizeStr = sprintf('%ux', sz);
    
    % Remove the last 'x' that is superfluous.
    sizeStr(end) = [];
    
