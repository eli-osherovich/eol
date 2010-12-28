function [val, grad, hess] = calc_EDx(x, Ax, FuncAxStruct, funcX, AD, D, ProjFlag, complexVarsFlag)
% CALC_EDX - calculates E(Dx) = F(ADx) + f(Dx)

% do not perform projection on V by default.
if nargin < 7
    ProjFlag = false;
end

% assume real variables by default.
if nargin < 8
    complexVarsFlag = false;
end

switch nargout
    case 1
        val1 = calc_FAx_series(FuncAxStruct, Ax, AD, ProjFlag, complexVarsFlag);
        val2 = calc_Fx_series(funcX, x, D, ProjFlag, complexVarsFlag);
        val = val1+val2;
    case 2
        [val1, grad1] = calc_FAx_series(FuncAxStruct, Ax, AD, ProjFlag, complexVarsFlag);
        [val2, grad2] = calc_Fx_series(funcX, x, D, ProjFlag, complexVarsFlag);
        val = val1 + val2;
        grad = grad1 + grad2;
    case 3
        [val1, grad1, hess1] = calc_FAx_series(FuncAxStruct, Ax, AD, ProjFlag, complexVarsFlag);
        [val2, grad2, hess2] = calc_Fx_series(funcX, x, D, ProjFlag, complexVarsFlag);
        val = val1 + val2;
        grad = grad1 + grad2;
        hess = hess1 + hess2;
    otherwise
        error('Opps');
end
     
