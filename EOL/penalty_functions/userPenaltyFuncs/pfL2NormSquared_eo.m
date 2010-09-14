function [val, grad, hessV, hessVfunc] = nonnegative_penalty_quadratic(x,V)
% NONEGATIVE PENALTY QUADRATIC - compute penalty for negative values 
% val = 1/2 * x^2  for x<0 and x belongs to params.support
% val = sv_weight * 1/2*x^2 for x not in the support
% val = 0   otherwise


% Copyright 2008-2010 Eli Osherovich.

assert(isvector(x));
assert(isreal(x));

%sv_weight = 1e5;
violations = false(size(x));
%support_violations = false(size(x));
%if isfield(params, 'support') && ~isempty(params.support) && all(size(x) == size(params.support))
%        violations(params.support) = x(params.support) < 0;
%        support_violations(~params.support) = x(~params.support)~=0;
%else
    violations = x<0;
%end

val = sum(x(violations).^2);
%val = val + sv_weight*sum(x(support_violations).^2)/2;

if nargout > 1, % gradient requested
    grad = zeros(size(x));
    grad(violations) = 2*x(violations);
    %grad(support_violations) = sv_weight*x(support_violations);
    
    
    switch nargout
        case 3
            if nargin == 1 || isempty(V), % return Hessian
                hessV = 2*diag(violations(:)); % + sv_weight*diag(support_violations(:));
            else
                hessV = hessMult(V);
            end
        case 4
            hessV = [];
            hessVfunc = @hessMult;
    end
end

    function hessV = hessMult(V)
        % preallocate space
        hessDiag = 2*violations(:);% + sv_weight*support_violations(:);
        
        if isnumeric(V) % we have one vector
            assert(isvector(V));
            hessV = hessDiag.*V(:);
        elseif iscell(V) % multiple vectors
            hessV   = cell(size(V));
            for i = 1:numel(V)
                if isempty(V{i})
                    hessV{i} = [];
                else
                    assert(isvector(V{i}));
                    tmp =  hessDiag.*V{i};
                    hessV{i} = tmp(:);
                end
            end
        else
            error('Wrong V: it must be either a vector or cell');
        end
    end
end

