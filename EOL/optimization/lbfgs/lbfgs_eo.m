function [x, funcVal, Output] = lbfgs_eo(x0, FuncAxStruct, ...
    funcX, Options)
% LBFGS - perform unconstrained minimization with the L-BFGS algorithm.


% Copyright 2010 Eli Osherovich.



% Save the original shape of x0. While working internally with a column
% vector, we will reshape x to its original shape on exit.
xOrigSize = size(x0);


%% Use default options if Options is not provided or empty.
if nargin < 4
    Options = struct();
end

%% Set parameters (default value)
[   x,...               % Initial x: obtained from x0 by converting
    ...                 % to column vector and changing type to
    ...                 % complex/real if necessary
    maxIter,...         % Maximal number of iterations (200)
    complexVarsFlag,... % Indicator whether the variables X are complex (false)
    nPrev,...           % Number of previous steps/gradients to remember (10)
    useMex,...          % Shall we use MEX files? (true)
    tolX, ...           % Step size tolerance (1e-8)
    tolFun, ...         % Function value tolerance (1e-8)
    targetFuncVal, ...  % Objective function target value (-realmax);
    tolGrad,...         % Gradient norm tolerance (1e-8)
    display,...         % Progress report (true)
    saveDir...          % Directory to save current x ('' = do not save)
    ] = lbfgsGetOptions_eo(x0, Options);


% Check if the third output argument (Output) was requested. 
if nargout > 2
    outputRequested = true;
    Output.allFuncVals = NaN(maxIter+1, 1);
    Output.allGradNorms = NaN(maxIter+1, 1);
else
    outputRequested = false;
end

%% Variables initialization and memory allocation.
% Preallocate space for previous steps and grads.
prevGrads = cell(1, nPrev);
prevSteps = cell(1, nPrev);

% Calculate Ax for the initial point.
Ax = applyMapping(FuncAxStruct, x);

% Calculate initial function value and gradient.
[funcVal, grad] = calcObjFunc(x, Ax, FuncAxStruct, funcX, complexVarsFlag);
gradNorm = norm(grad);

% Initialize counters.
funcCount = 1;
iter = 0;

% Print initial state (if requested).
if display
	fprintf('%10s %10s %15s %15s %15s\n','Iteration','FunEvals','Step Length','Function Val','Opt Cond');
    fprintf('%10d %10d %15.5e %15.5e %15.5e\n', iter, ...
        funcCount, 0, funcVal, gradNorm);
end
% Save current x (if requested).
if ~isempty(saveDir)
    save(fullfile(saveDir, int2str(iter)), 'x');
end

% Test termination criteria.
[done, exitFlag, exitMsg] = testTermCriteria(...
    iter, x, grad, funcVal, ...
    Inf, gradNorm, Inf,...
    [], maxIter, tolX, tolGrad, tolFun, targetFuncVal, complexVarsFlag);

% If Output structure was requested update the list of all function
% values and gradient norms (used to review/visualize minimization
% progress).
if outputRequested 
    Output.allFuncVals(iter + 1) = funcVal;
    Output.allGradNorms(iter + 1) = gradNorm;
end


% Initial step and search direction scale factor.
t0 =  1;
H0 = min(1, 1/gradNorm);


validIdx = 0;
wrapAround = 0;

%% Run iterations.
while ~done
    % Update iteration counter.
    iter = iter + 1;
    
    % Calculate new search direction.
    if useMex
        d = lbfgsDirPersistC_eo(validIdx, wrapAround, -grad, ...
            prevSteps, prevGrads, H0);
    else
        d = lbfgsDirPersist_eo(validIdx, wrapAround, -grad, ...
            prevSteps, prevGrads, H0);
    end
    
    % Calculate d's linear transform (used by line search).
    %Ad = applyMapping(FuncAxStruct, d);
    Ad = cell(size(Ax));
    
    % Save current point data.
    gradOld = grad;
    funcValOld = funcVal;
    
    % Peform line search.
    [t, x, funcVal, grad, Ax, LSoutput] =  wolfeLS_eo(t0, x, funcVal, grad, ...
        d, Ax, Ad, FuncAxStruct, funcX, Options);
    funcCount = funcCount + LSoutput.funcCount;
    
    
    % Go to the new point.
    gradNorm = norm(grad);
    step = t*d;
    stepNorm = norm(step);
    
    
    
    % Update LBFGS memory.
    [prevSteps, prevGrads, H0, validIdx, wrapAround] = ...
        lbfgsUpdate_eo(validIdx, wrapAround, step, grad-gradOld,...
        prevSteps, prevGrads, H0);

    % Print progress (if requested).
    if display
        fprintf('%10d %10d %15.5e %15.5e %15.5e\n', iter, ...
            funcCount, t, funcVal, gradNorm);
    end
    % If Output structure was requested update the list of all function
    % values and gradient norms (used to review/visualize minimization
    % progress).
    if outputRequested
        Output.allFuncVals(iter + 1) = funcVal;
        Output.allGradNorms(iter + 1) = gradNorm;
    end
    % Save current x (if requested).
    if ~isempty(saveDir)
        save(fullfile(saveDir, int2str(iter)), 'x');
    end
    
    % Test termination criteria.
    [done, exitFlag, exitMsg] = testTermCriteria(...
        iter, x, grad, funcVal, ...
        stepNorm, gradNorm, abs(funcVal - funcValOld),...
        LSoutput, maxIter, tolX, tolGrad, tolFun, targetFuncVal, complexVarsFlag);
   
end
% Display exit message (if requested).
if display
    fprintf('%s\n', exitMsg);
end

% Reshape x to its original size.
x = reshape(x, xOrigSize);

% Update Output structure (if requested).
if outputRequested
    % Set output struct fields.
    Output.finalGrad = grad;
    Output.nIterations = iter;
    Output.funcCount = funcCount;
    Output.exitMsg = exitMsg;
    Output.exitFlag = exitFlag;
    
    % Clip unused function values/gradient norms.
    Output.allFuncVals(iter + 2:end) = [];
    Output.allGradNorms(iter + 2:end) = [];
    
    % Make sure that we use LSoutput only if the linesearch was called at least
    % once.
    if iter > 0
        Output.lsExitMsg = LSoutput.exitMsg;
        Output.lsExitFlag = LSoutput.exitFlag;
    else
        Output.lsExitMsg = '';
        Output.lsExitFlag = NaN;
    end
end


function [done, exitFlag, exitMsg] = testTermCriteria(...
    iter, x, grad, funcVal, ...
    stepNorm, gradNorm, funcValDiff,...
    LSoutput, maxIter, tolX, tolGrad, tolFun, targetFuncVal, complexVarsFlag)

done = false;
exitFlag = 0;
exitMsg = '';

% Check if x has valid values.
if ~isValid(x, complexVarsFlag)
    done = true;
    exitMsg = 'Invalid X values';
    exitFlag = -1;
    return;
end

% Check if grad has valid values.
if ~isValid(grad, complexVarsFlag)
    done = true;
    exitMsg = 'Invalid Gradient values';
    exitFlag = -2;
    return;
end

% Check if fVal is valid.
if ~isValid(funcVal)
    done = true;
    exitMsg = 'Invalid Function value';
    exitFlag = -3;
    return;
end

% Check whether line-search routine finished successfully. This test
% skipped if we are checking initial conditions, i.e., iter == 0.
if iter > 0 
    if LSoutput.exitFlag < 0
        % An error occured during line-search.
        done = true;
        exitFlag = -4;
        exitMsg = sprintf('Unrecoverable error in the line-search routine:\n %s\n', ...
            LSoutput.exitMsg);
        return;
    elseif LSoutput.exitFlag > 0
        % Line-search failed to find a suitable step-length. Maybe required
        % tolerance is too strict.
        done = true;
        exitFlag = -4;
        exitMsg = sprintf('Line-search cound not find a suitable step-length:\n %s\n', ...
            LSoutput.exitMsg);
        return;
    end
end
    

% If gradient's norm is small (below tolGrad) we probably found a local
% minimum. Hence, exit status is zero which indicates a success.
if gradNorm <= tolGrad
    done = true;
    exitFlag = 0;
    exitMsg = sprintf(['Gradient magnitude is smaller than tolGrad (%g) ' ...
                       'tolerance.'], tolGrad);
    return;
end

% If we reached the maximum allowed number of iterations terminate with
% exit status set to one.
if iter >= maxIter
    done = true;
    exitFlag = 1;
    exitMsg = sprintf(['Number of iterations (%d) exceeded maximum ' ...
                       'allowed.'], maxIter);
    return;
end

if stepNorm <= tolX
    done = true;
    exitFlag = 2;
    exitMsg = sprintf('Step length is shorter than tolX (%g) tolerance.', ...
                  tolX);
    return;
end

if funcValDiff <= tolFun
    done = true;
    exitFlag = 3;
    exitMsg = sprintf(['Change in the objective function value is smaller ' ...
                       'than tolFun (%g) tolerance.'], tolFun);
    return;
end

if funcVal <= targetFuncVal
    done = true;
    exitFlag = 4;
    exitMsg = sprintf(['Objective function value has reached the desired ' ...
                       'target targetFuncVal (%g).'], targetFuncVal);
end
