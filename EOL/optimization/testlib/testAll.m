function testAll(dirName, fileOutName)
%TESTALL - test all minimization problems in a directory.
% 
% Usage:
%-------
% TESTALL(DIRNAME) - test all minimization problems in the directory
% DIRNAME. Each problem is implemented in one shared library (.so or .dll
% file).
% The function tests correctness of Jacobian/Gradient calculations.
%
% IMPLEMENTATION DETAILS:
%------------------------
% See objFunAll()



% Copyright 2010-2011 Eli Osherovich.

    
    
% Use current directory if not specified.
if 0 == nargin
    dirName = pwd;
end

% Use standard output if fileOutName was not provided.
if nargin < 2 || isempty(fileOutName)
    fileID = 1;
else
    fileID = fopen(fileOutName, 'w');
end
   
% All shared libraries in the directory.
allSlibs = dir(fullfile(dirName, '*.so'));


% Print header.
fprintf(fileID, '%-15s %-15s %-15s %-15s %-15s\n', 'Function', ['Jac Abs ' ...
                    'Err'], 'Jac Rel Err', 'Grad Abs Err', 'Grad Rel Err');
fprintf(fileID, ['-------------------------------------------------------' ...
                 '------------------------\n']); 
for i = 1:numel(allSlibs)
    libFile = fullfile(dirName, allSlibs(i).name);
    [~, funcName] = fileparts(libFile);
    
    % Get the standard starting point (defined by the problem).
    [~, ~, ~, ~, x] = objFuncAll(libFile);
    
    % Generate random point near the standard one.
    x = x + rand(size(x));
    
    % Calculate analytical Jacobian and gradient.
    [~, ~, fj, g] = objFuncAll(libFile, x);
    
        
    % Calculate numerical Jacobian.
    N_Jac = calcNumJacobian_eo(x, @objFuncSOE, 'precise');
    
    % Calculate numerical gradient.
    N_Grad = (calcNumJacobian_eo(x, @objFuncSOS, 'precise'))';
    
    % Jacobian max errors.
    jacAbsErr = max(abs(fj(:) - N_Jac(:)));
    jacRelErr = max(abs((fj(:) - N_Jac(:))./fj(:)));
    
    % Gradient max errors.
    gradAbsErr = max(abs(g(:) - N_Grad(:)));
    gradRelErr = max(abs((g(:) - N_Grad(:))./g(:)));
    
    % Prevent NaNs in relative errors when abs error is zero.
    if 0 == jacAbsErr 
        jacRelErr = 0;
    end
    if 0 == gradAbsErr
        gradRelErr = 0;
    end
    
    % Print results.
    fprintf(fileID, '%-15s %-15e %-15e %-15e %-15e\n', funcName, ...
        jacAbsErr, jacRelErr, gradAbsErr, gradRelErr);
end
% Dummy system of equations problem
% for numerical Jacobian calculation.
    function fvec = objFuncSOE (x)
        fvec = objFuncAll(libFile, x);
    end

% Dummy sum-of-squares problem
% for numerical gradient calculation.
    function fval = objFuncSOS (x)
        [~, fval] = objFuncAll(libFile, x);
    end
end
