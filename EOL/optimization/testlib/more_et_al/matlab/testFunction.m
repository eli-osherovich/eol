function [absErr, relErr] = testFunction(func)

% get standard starting point 
[~, ~, ~, ~, x] = func();

x = x + rand(size(x));

% compute analytical Jacobian
[~, A_Jac] = func(x);

% calculate numerical Jacobian
N_Jac = (calcNJacCDExt_eo(func, x))';

% absolute and relative errors
absErr = max(abs(A_Jac(:) - N_Jac(:)));
relErr = max(abs((A_Jac(:) - N_Jac(:))./A_Jac(:)));

if relErr > 1e-8
    disp(A_Jac-N_Jac)
end
