function benchUO (dirName, fileOutName)
%BENCHUO - unconstrained optimization benchmarks

    
% Copyright 2010 Eli Osherovich.
    

% Use current directory by default.
if 0 == nargin || isempty(dirName)
    dirName = pwd;
end

% Use standard output if fileOutName was not provided.
if nargin < 2 || isempty(fileOutName)
    fileID = 1;
else
    fileID = fopen(fileOutName, 'w');
end

% all shared libraries in the directory
allSlibs = dir(fullfile(dirName, '*.so'));


% print header
fprintf(fileID, '%-10s %-15s %-15s %-7s %-7s %-9s %-s\n', 'Function', 'F Value', ...
        'G Norm', 'Iter', 'Fcnt', 'Exit', 'Algorithm');
fprintf(fileID, ['---------------------------------------------------------------' ...
         '----------------\n']); 

     
%% generate common optimization parameters 
options = optimset();
options = optimset(options, 'Jacobian', 'on');
options = optimset(options, 'GradObj', 'on');
options = optimset(options, 'MaxIter', 1000);
options = optimset(options, 'TolX', 1e-10);
options = optimset(options, 'TolFun', 1e-10);
options = optimset(options, 'Display', 'off');

     
% run over all problems (shared libraries)
for i = 1:numel(allSlibs)
    fprintf(fileID, '\n');
    
    libFile = fullfile(dirName, allSlibs(i).name);
    [~, funcName] = fileparts(libFile);
    
    % get standard starting point 
    [~, ~, ~, ~, x0] = objFuncAll(libFile);
    
%     %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % MATLAB's fuminunc: Medium-Scale mode              %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     options = optimset(options, 'LargeScale', 'off');
%     
%     [~,fval,exitflag,output] = fminunc(@objFuncSOS, x0, options);
%     % print results
%     fprintf(fileID, '%-10s %-15g %-15g %-7d %-7d %-9s %-s\n', funcName, ...
%             fval, output.firstorderopt, output.iterations, output.funcCount, ...
%             status2str(exitflag), output.algorithm );
%     
%     %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %   MATLAB's fuminunc: Large-Scale mode             %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     options = optimset(options, 'LargeScale', 'on');
%         
%     [~, fval, exitflag, output] = fminunc(@objFuncSOS, x0, options);
%     % print results
%     fprintf(fileID, '%-10s %-15g %-15g %-7d %-7d %-9s %-s\n', funcName, ...
%             fval, output.firstorderopt, output.iterations, output.funcCount, ...
%             status2str(exitflag), output.algorithm );
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          MINFUNC                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %options.Display = 'on';
    [~, fval, exitflag, output] = minFunc(@objFuncSOS, x0, options);
    fprintf(fileID, '%-10s %-15g %-15g %-7d %-7d %-9s %-s\n', funcName, ...
            fval, output.firstorderopt, output.iterations, output.funcCount, ...
            status2str(exitflag), 'minFunc' );
        
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         LBFGS_EO                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fxStruct = struct('function', @objFuncSOS);
    %options.Display = 'on';
    [~, fval, exitflag, output] = lbfgsO_eo(x0, [], fxStruct, options);
    fprintf(fileID, '%-10s %-15g %-15g %-7d %-7d %-9s %-s\n', funcName, ...
            fval, output.firstorderopt, output.iterations, output.funcCount, ...
            status2str(exitflag), 'lbfgsO_eo' );
        
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         LBFGS_EO  - TEST                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fxStruct = struct('function', @objFuncSOS);
    [~, fval, output] = lbfgs_eo(x0, [], fxStruct, options);
    fprintf(fileID, '%-10s %-15g %-15g %-7d %-7d %-9s %-s\n', funcName, ...
            fval, output.firstOrderOpt, output.nIterations, output.funcCount, ...
            status2str(output.exitFlag), 'lbfgs_eo' );
        %break;
%     %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % non-linear least squares (trust region reflective)  %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     options = optimset(options, 'Algorithm', 'trust-region-reflective' );
%     
%     [~, resnorm, ~, exitflag, output] = lsqnonlin(@objFuncNLLS, x0, [], ...
%                                                   [], options);
%     fprintf(fileID, '%-10s %-15g %-15g %-7d %-7d %-9s %-s\n', funcName, ...
%             resnorm, output.firstorderopt, output.iterations, output.funcCount, ...
%             status2str(exitflag), output.algorithm );
%     
%     %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % non-linear least squares (Levenberg-Marquardt)      %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     options = optimset(options, 'Algorithm', 'levenberg-marquardt');
%     
%     [~, resnorm, ~, exitflag, output] = lsqnonlin(@objFuncNLLS, x0, [], ...
%                                                   [], options);
%     fprintf(fileID, '%-10s %-15g %-15g %-7d %-7d %-9s %-s\n', funcName, ...
%             resnorm, output.firstorderopt, output.iterations, output.funcCount, ...
%             status2str(exitflag), output.algorithm );

end

% Close file if it is not the standard output
if 1 ~= fileID
    fclose(fileID);
end

% dummy sum-of-squares problem
function [fval, grad]  = objFuncSOS (x, ~)
    [~, fval, ~, grad] = objFuncAll(libFile, x);
end

% dummy non-linear least squares problem
function [fvec, jac] = objFuncNLLS(x)
    [fvec, ~, jac] = objFuncAll(libFile, x);
end
end


function str = status2str(status)
    switch status
      case -4
          str = 'LS Error';
          
        case -3
            str = 'F Error';
            
        case -2
            str = 'G Error';
            
        case -1
            str = 'X Error';
            
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
        str = '???';     % Unknown.
    end
end
