function [val, grad] = minFunc_wrapper(x, func_Ax_Struct, func_x_Struct, ComplexVarsFlag)

% assume real variables by default
if nargin < 4
    ComplexVarsFlag = false;
end


if ComplexVarsFlag
    n = numel(x);
    x = reshape(complex(x(1:n/2), x(n/2+1:end)), n/2, 1);
end

% generate empty cell array of proper size (used by some functions)
empty = cell(size(func_Ax_Struct));
   
Ax = calculate_linop('forward', func_Ax_Struct, x);

switch nargout
     case 1
         val = calc_EDx(x, Ax, func_Ax_Struct, func_x_Struct, empty, [], false, ComplexVarsFlag);
     case 2
         [val, grad] = calc_EDx(x, Ax, func_Ax_Struct, func_x_Struct, empty, [], false, ComplexVarsFlag);
         
         if ComplexVarsFlag
             grad = [real(grad); imag(grad)];
         end
 end
    
       
