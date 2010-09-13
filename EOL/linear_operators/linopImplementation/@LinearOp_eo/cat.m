function newObj = cat(dim, varargin)



% Copyright 2010 Eli Osherovich.



% Check that the concatenation dimenstion is valid.
validateattributes(dim, {'numeric'}, {'integer', ...
    'positive', 'scalar'});

% Currently only two-dimensional arrays (matices) are allowed.
% Higher dimenstional array (tensors) are not supported at the moment.
switch dim
    case 1
        newObj = LinearOpMatrix_eo(varargin');
    case 2
        newObj = LinearOpMatrix_eo(varargin);
    otherwise
            error('EOL:LinearOp:cat:WrongArg', ...
                'Only two dimensional arrays are supported.');
end
        
