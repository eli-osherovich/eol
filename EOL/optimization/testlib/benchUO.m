function benchUO (dirName)
%BENCHUO - unconstrained optimization benchmarks

    
% Copyright 2010 Eli Osherovich.
    

% use current directory by default
if 0 == nargin || isempty(dirName)
    dirName = pwd;
end

% all shared libraries in the directory
allSlibs = dir(fullfile(dirName, '*.so'));


% print header
fprintf('%-10s %-15s %-15s %-7s %-7s %-9s %-s\n', 'Function', 'F Value', ...
        'G Norm', 'Iter', 'Fcnt', 'Exit', 'Algorithm');
fprintf(['---------------------------------------------------------------' ...
         '----------------\n']); 

% run over all problems (shared libraries)
for i = 1:numel(allSlibs)
    fprintf('\n');
    
    libFile = fullfile(dirName, allSlibs(i).name);
    [~, funcName] = fileparts(libFile);
    
    % get standard starting point 
    [~, ~, ~, ~, x0] = objFuncAll(libFile);
    
    %% MATLAB's fuminunc: Medium-Scale mode
    options = optimset();
    options = optimset(options, 'LargeScale', 'off');
    options = optimset(options, 'GradObj', 'on');
    options = optimset(options, 'Display', 'off');
    
    [~,fval,exitflag,output] = fminunc(@objFuncSOS, x0, options);
    % print results
    fprintf('%-10s %-15g %-15g %-7d %-7d %-9s %-s\n', funcName, ...
            fval, output.firstorderopt, output.iterations, output.funcCount, ...
            status2str(exitflag), output.algorithm );
    
    %% MATLAB's fuminunc: Large-Scale mode
    options = optimset();
    options = optimset(options, 'LargeScale', 'on');
    options = optimset(options, 'GradObj', 'on');
    options = optimset(options, 'Display', 'off');
    
    [~, fval, exitflag, output] = fminunc(@objFuncSOS, x0, options);
    % print results
    fprintf('%-10s %-15g %-15g %-7d %-7d %-9s %-s\n', funcName, ...
            fval, output.firstorderopt, output.iterations, output.funcCount, ...
            status2str(exitflag), output.algorithm );
    
    %% minFunc
    options = optimset();
    options = optimset(options, 'Display', 'off');
    [~, fval, exitflag, output] = minFunc(@objFuncSOS, x0, options);
    fprintf('%-10s %-15g %-15g %-7d %-7d %-9s %-s\n', funcName, ...
            fval, output.firstorderopt, output.iterations, output.funcCount, ...
            status2str(exitflag), 'minFunc' );
    
    %% non-linear least squares (trust region reflective)
    options = optimset('Jacobian', 'on');
    options = optimset(options, 'Display', 'off');
    [~, resnorm, ~, exitflag, output] = lsqnonlin(@objFuncNLLS, x0, [], ...
                                                  [], options);
    fprintf('%-10s %-15g %-15g %-7d %-7d %-9s %-s\n', funcName, ...
            resnorm, output.firstorderopt, output.iterations, output.funcCount, ...
            status2str(exitflag), output.algorithm );
    
    %% non-linear least squares (Levenberg-Marquardt)
    options = optimset('Jacobian', 'on');
    options = optimset(options, 'Display', 'off');
    options = optimset(options, 'Algorithm', 'levenberg-marquardt');
    [~, resnorm, ~, exitflag, output] = lsqnonlin(@objFuncNLLS, x0, [], ...
                                                  [], options);
    fprintf('%-10s %-15g %-15g %-7d %-7d %-9s %-s\n', funcName, ...
            resnorm, output.firstorderopt, output.iterations, output.funcCount, ...
            status2str(exitflag), output.algorithm );

end

% dummy sum-of-squares problem
function [fval, grad]  = objFuncSOS (x)
    [~, fval, ~, grad] = objFuncAll(libFile, x);
end

% dummy non-linear least squares problem
function [fvec, jac] = objFuncNLLS(x)
    [fvec, ~, jac] = objFuncAll(libFile, x);
end
end


function str = status2str(status)
    switch status
      case 0 
        str = 'MaxIter';  % maximum number of iterations reached
                          
      case 1 
        str = 'SmallGrd'; % magnitude of gradient is smaller than tolFun tolerance
                          % converged in case of lsqnonlin  
        
      case 2 
        str = 'SmallStp'; % change in x was smaller than tolX tolerance
        
      case 3 
        str = 'SmallDcr'; % change in the objective function (or
                          % residual) was less than specified tolerance
        
      case 4 % lsqnonlin only
        str = 'SmallDD';  % small directional derivative
        
      case 5 
        str = 'PredDcr';  % Predicted decrease in the objective function was less
                          % than the TolFun tolerance
        
      otherwise 
        str = '???';     % unkonown
    end
end
