function [val, grad] = minFunc_wrapper(x, FuncAxStruct, funcX, complexVarsFlag)

% Assume real variables by default.
if nargin < 4
    complexVarsFlag = false;
end


if complexVarsFlag
    n = numel(x);
    x = reshape(complex(x(1:n/2), x(n/2+1:end)), n/2, 1);
end

 
% Calculate Ax for the initial point.
Ax = applyMapping(FuncAxStruct, x);

switch nargout
     case 1
         
         val = calcObjFunc(x, Ax, FuncAxStruct, funcX, complexVarsFlag);
     
    case 2

         [val, grad] = calcObjFunc(x, Ax, FuncAxStruct, funcX, complexVarsFlag);
         if complexVarsFlag
             grad = [real(grad); imag(grad)];
         end
end
    
       
