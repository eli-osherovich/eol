function [val, grad, hess] = calc_FAx_series(fAxStruct, Ax, Av, ProjFlag, ComplexVarsFlag)


% (c) Copyright 2008-2010 Eli Osherovich.


% fast exit if empty problem
if isempty(fAxStruct)
    val = 0;
    grad = 0;
    hess = 0;
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
        for i = 1:length(fAxStruct)
       
            weight = find_weight(fAxStruct(i));

            val = val + weight * calc_Fx_series(fAxStruct(i).fxStruct, Ax{i}, [], ProjFlag, ComplexVarsFlag);
        end


    case 2
        val = 0; grad = 0;
        
        for i = 1:length(fAxStruct)
            weight = find_weight(fAxStruct(i));
            
            % Attention!
            % this is probably not what you expect!
            % we calculate dot product between the gradient and v
            % if v was given (non-empty)
            [val_tmp, grad_tmp] = calc_Fx_series(fAxStruct(i).fxStruct, Ax{i}, Av{i}, ProjFlag, ComplexVarsFlag);
        
            val = val + weight * val_tmp;
            
            if ~ProjFlag || isempty(Av{i})
                grad = grad + weight * calculate_linop('adjoint', fAxStruct, grad_tmp, i);
            else
                grad = grad + weight * grad_tmp;
            end
        end
        
        % make sure the gradient is real if opt. variables are real
        if ~ComplexVarsFlag
            grad = real(grad);
        end
        
    case 3
        val = 0; 
        grad = 0; 
        hess = 0;
        
        for i = 1:length(fAxStruct)
            
            weight = find_weight(fAxStruct(i));
            
            [val_tmp, grad_tmp, hessAv_tmp] = calc_Fx_series(fAxStruct(i).fxStruct, Ax{i}, Av{i}, ProjFlag, ComplexVarsFlag);
            
            val = val + weight * val_tmp;
            
            if ~ProjFlag || isempty(Av{i})
                grad = grad + weight * calculate_linop('adjoint', fAxStruct, grad_tmp, i);
                % hessian is calculated by formula A'*H*A = A'*H*A*I
                % it may be extremely slow.
                hess = hess + weight * (...
                    calculate_linop('adjoint', fAxStruct, hessAv_tmp, i) * ...
                    calculate_linop('forward', fAxStruct, eye(numel(Ax{i})), i));
            else
                grad = grad + weight * grad_tmp;
                hess = hess + weight * hessAv_tmp;
            end
        end
        
        % make sure the gradient and Hessian are real if opt. variables are real
        if ~ComplexVarsFlag
            grad = real(grad);
            hess = real(hess);
        end
        
    otherwise
        error('Function returns one to three outputs');
end





