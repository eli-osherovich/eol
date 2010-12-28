function [val, grad] = minFunc_wrapper(x, FuncAxStruct, funcX, complexVarsFlag)

% assume real variables by default
if nargin < 4
    complexVarsFlag = false;
end


if complexVarsFlag
    n = numel(x);
    x = reshape(complex(x(1:n/2), x(n/2+1:end)), n/2, 1);
end

% generate empty cell array of proper size (used by some functions)
empty = cell(size(FuncAxStruct));
   
Ax = calculate_linop('forward', FuncAxStruct, x);

switch nargout
     case 1
         val = calc_EDx(x, Ax, FuncAxStruct, funcX, empty, [], false, complexVarsFlag);
     case 2
         [val, grad] = calc_EDx(x, Ax, FuncAxStruct, funcX, empty, [], false, complexVarsFlag);
         
         if complexVarsFlag
             grad = [real(grad); imag(grad)];
         end
 end
    
       
