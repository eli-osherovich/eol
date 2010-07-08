function testAll(dirName)
%TESTALL - test all minimization problems in a directory.
% 
% Usage:
%-------
% TESTALL(DIRNAME) - test all minimization problems in the directory
% DIRNAME. Each problem is implemented in one shared library (.so or .dll
% file).
%
% IMPLEMENTATION DETAILS:
%------------------------
% See objFunAll()
    
    
% Copyright 2010 Eli Osherovich.

% use current directory if not sepcified
if 0 == nargin
    dirName = pwd;
end
    
% all shared libraries in the directory
allSlibs = dir(fullfile(dirName, '*.so'));


% print header
fprintf('%15s %15s %15s %15s %15s\n', 'Function',...
    'Jac Abs Err', 'Jac Rel Err', ...
    'Grad Abs Err', 'Grad Rel Err');

for i = 1:numel(allSlibs)
    libFile = fullfile(dirName, allSlibs(i).name);
    [~, funcName] = fileparts(libFile);
    
    % get standard starting point 
    [~, ~, ~, ~, x] = objFuncAll(libFile);
    
    % generate random point near the standard one
    x = x + rand(size(x));
    
    % calculate analytical Jacobian and grient
    [~, ~, fj, g] = objFuncAll(libFile, x);
    
        
    % calculate numerical Jacobian
    N_Jac = (calcNJacCDExt_eo(@objFuncSOE, x))';
    
    % calculate numerical gradient
    N_Grad = (calcNJacCDExt_eo(@objFuncSOS, x))';
    
    % Jacobian max errors
    jacAbsErr = max(abs(fj(:) - N_Jac(:)));
    jacRelErr = max(abs((fj(:) - N_Jac(:))./fj(:)));
    
    % gradient max errors
    gradAbsErr = max(abs(g(:) - N_Grad(:)));
    gradRelErr = max(abs((g(:) - N_Grad(:))./g(:)));
    
    % prevent NaNs in relative error
    if 0 == jacAbsErr 
        jacRelErr = 0;
    end
    if 0 == gradAbsErr
        gradRelErr = 0;
    end
    
    % print results
    fprintf('%15s %15g %15g %15g %15g\n', funcName, ...
        jacAbsErr, jacRelErr, gradAbsErr, gradRelErr);
end
% dummy system of equations problem
% for numerical Jacobian calculation
    function fvec = objFuncSOE (x)
        fvec = objFuncAll(libFile, x);
    end

% dummy sum-of-squares problem
% for numerical gradient calculation
    function fval = objFuncSOS (x)
        [~, fval] = objFuncAll(libFile, x);
    end
end
