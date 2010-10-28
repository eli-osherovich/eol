function [val, grad, hessV] = subsref(self, restArgs)

% Make sure that the function was called with the regular parentheses: ()

%if ~strcmp(restArgs(1).type, '()')
if '(' ~= restArgs(1).type(1)
    val = builtin('subsref', self, restArgs);
    return;
end

x = restArgs.subs{1};

% By the convention x is assumed to be a column vector.
validateattributes(x, {'numeric'}, {'column'});

mFactor = self.multFactor;

switch nargout

    case {0, 1}
        % Return the value at x.
        val = doCalculations(self, x);
         
    case 2
        % Return the value and gradient.
        [val, grad] = doCalculations(self, x);
         
    case 3
        v = restArgs.subs{2};
        % Return the value, gradient and Hessian times vector.
        [val, grad, hessMultVectorFunc] = doCalculations(self, x);
         
        if isnumeric(v) % one vector
            assert(isvector(v));
            hessV = hessMultVectorFunc(v);
             
        elseif iscell(v) % multiple vectors
            % preallocate space
            hessV = cell(size(v));
            for i = 1:numel(v)
                if ~isempty(v{i})
                    assert(isvector(v{i}));
                    hessV{i} = hessMultVectorFunc(v{i});
                 end
            end
        else
            error('EOL:PenaltyFunc:WrongArgType', ...
                'V must be a numeric vector of cell array');
        end
        
    otherwise
        error('EOL:PenaltyFunc:WrongOutNum', ...
            'The function retrurns up to 3 outputs');
end

% Apply mFactor (if not unity).
if mFactor ~= 1
    val = val * mFactor;
    grad = grad * mFactor;
    if isnumeric(hessV)
        hessV = hessV * mFactor;
    else % hessV is a cell
        for i = 1:numel(hessV)
            hessV{i} = hessV{i} * mFactor;
        end
    end
end
