function [val, gradV, VhessV] = calc_Fx_series(fxStruct, x, V, ProjFlag, ComplexVarsFlag)


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
    ComplexVarsFlag = false;
end

switch nargout
    case 1
        val = 0;
        for i = 1:length(fxStruct)
            
            weight = find_weight(fxStruct(i));
            
            if isfield(fxStruct(i), 'func_args')
                val = val + weight * fxStruct(i).function(x, [], fxStruct(i).func_args{:});
            else
                val = val + weight * fxStruct(i).function(x, []);
            end
        end
        
    case 2
        val = 0;
        grad = zeros(size(x));
        
        for i = 1:length(fxStruct)
                        
            weight = find_weight(fxStruct(i));

            if isfield(fxStruct(i), 'func_args')
                [val_tmp, grad_tmp] = fxStruct(i).function(x, [], fxStruct(i).func_args{:});
            else
                [val_tmp, grad_tmp] = fxStruct(i).function(x, []);
            end
            val = val + weight * val_tmp;
            grad = grad + weight * grad_tmp;
        end
        
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
            if ~ComplexVarsFlag
                gradV = real(gradV);
            end
        else
            gradV = grad;
        end
	
	            
        
    case 3
       val = 0; 
       grad = zeros(numel(x), 1);
       hessV = cell(size(V));
       [hessV{:}] = deal(zeros(numel(x),1));
        
        for i = 1:length(fxStruct)
                                    
            weight = find_weight(fxStruct(i));

            if isfield(fxStruct(i), 'func_args')
                [val_tmp, grad_tmp, hessV_tmp] = fxStruct(i).function(x, V, fxStruct(i).func_args{:});
            else
                [val_tmp, grad_tmp, hessV_tmp] = fxStruct(i).function(x,V);
            end
            val = val + weight * val_tmp;
            grad = grad + weight * grad_tmp;
            for k = 1:numel(hessV)
                if ~isempty(hessV_tmp{k})
                    hessV{k} = hessV{k} + weight * hessV_tmp{k};
                end
            end
        end
        
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
            if ~ComplexVarsFlag
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
