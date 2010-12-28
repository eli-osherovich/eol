function [val, gradV, VhessV] = calc_Fx_series(fxStruct, x, V, ProjFlag, complexVarsFlag)


% (c) Copyright 2008-2010 Eli Osherovich.


% fast exit if empty problem
if isempty(fxStruct)
    val = 0;
    gradV = 0;
    VhessV = 0;
    return
end

% do not perform projection on V by default.
if nargin < 4
    ProjFlag = false;
end

% assume real variables by default.
if nargin < 5
    complexVarsFlag = false;
end

switch nargout
    case 1
        val = fxStruct(x);

    case 2
       [val, grad] = fxStruct(x);
        
        % shall we calculate grad'*V ?
        if ProjFlag && ~isempty(V)
            if isnumeric(V)
                assert(isvector(V));
                gradV = V'*grad;
            elseif iscell(V)
                non_empty_idx = find(~cellfun('isempty',V));
                nnon_empty = sum(~cellfun('isempty',V));
                gradV = zeros(nnon_empty,1);
                for i = 1:nnon_empty
                    gradV(i) = V{non_empty_idx(i)}'*grad;
                end
            else
                error('Wrong V: it must be either a vector or cell');
            end
            % make sure the gradient is real if opt. variables are real
            if ~complexVarsFlag
                gradV = real(gradV);
            end
        else
            gradV = grad;
        end
	
	            
        
    case 3
       [val, grad, hessV] = fxStruct(x, V);
        
        if ProjFlag && ~isempty(V)
            non_empty_idx = find(~cellfun('isempty',V));
            nnon_empty = sum(~cellfun('isempty',V));
            gradV = zeros(nnon_empty,1);
            VhessV = zeros(nnon_empty);
            for i = 1:nnon_empty
                gradV(i) = V{non_empty_idx(i)}'*grad;
                for j=i:nnon_empty
                    VhessV(i,j) = V{non_empty_idx(i)}'*hessV{non_empty_idx(j)};
                    VhessV(j,i) = VhessV(i,j);
                end
            end
            % make sure the gradient and Hessian are real if opt. variables are real
            if ~complexVarsFlag
                gradV = real(gradV);
                VhessV = real(VhessV);
            end
        else
            gradV = grad;
            VhessV = hessV;
        end
	
        
    otherwise
        error('Function returns one to three outputs');
end
