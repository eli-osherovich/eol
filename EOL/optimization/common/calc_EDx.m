function [val, grad, hess] = calc_EDx(x, Ax, func_Ax_Struct, func_x_Struct, AD, D, ProjFlag, ComplexVarsFlag)
% CALC_EDX - calculates E(Dx) = F(ADx) + f(Dx)

% do not perform projection on V by default.
if nargin < 7
    ProjFlag = false;
end

% assume real variables by default.
if nargin < 8
    ComplexVarsFlag = false;
end

switch nargout
    case 1
        val1 = calc_FAx_series(func_Ax_Struct, Ax, AD, ProjFlag, ComplexVarsFlag);
        val2 = calc_Fx_series(func_x_Struct, x, D, ProjFlag, ComplexVarsFlag);
        val = val1+val2;
    case 2
        [val1, grad1] = calc_FAx_series(func_Ax_Struct, Ax, AD, ProjFlag, ComplexVarsFlag);
        [val2, grad2] = calc_Fx_series(func_x_Struct, x, D, ProjFlag, ComplexVarsFlag);
        val = val1 + val2;
        grad = grad1 + grad2;
    case 3
        [val1, grad1, hess1] = calc_FAx_series(func_Ax_Struct, Ax, AD, ProjFlag, ComplexVarsFlag);
        [val2, grad2, hess2] = calc_Fx_series(func_x_Struct, x, D, ProjFlag, ComplexVarsFlag);
        val = val1 + val2;
        grad = grad1 + grad2;
        hess = hess1 + hess2;
    otherwise
        error('Opps');
end
     
