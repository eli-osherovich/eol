function [val, grad] = calcObjFunc(x, Ax, FuncAxStruct, funcX, complexVarsFlag)
    % CALCOBJFUNC- calculates objective function F(Ax) + f(x)
    
    
    
    % Copyright 2010 Eli Osherovich.
    
    
    
    % Assume real variables by default.
    if nargin < 5
        complexVarsFlag = false;
    end
    
    % Preset outputs to zeros.
    val = 0;
    grad = 0;
    
    % Calculate the contribution of the first term.
    switch nargout
        
        case {0,1}
            if ~isempty(funcX)
                val = funcX(x);
            end
            
            for i = 1:length(FuncAxStruct)
                val_tmp = FuncAxStruct(i).function(Ax{i});
                val = val + val_tmp;
            end
            
        case 2
            if ~isempty(funcX)
                [val, grad] = funcX(x);
            end
            
            for i = 1:length(FuncAxStruct)
                [val_tmp, grad_tmp] = FuncAxStruct(i).function(Ax{i});
                
                val = val + val_tmp;
                
                grad_tmp = FuncAxStruct(i).linop' * grad_tmp;
                % make sure the gradient is real if opt. variables are real
                if ~complexVarsFlag
                    grad_tmp = real(grad_tmp);
                end
                grad = grad + grad_tmp;
            end
    end
