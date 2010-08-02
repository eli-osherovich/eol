function [x, funcVal, Output] = lbfgs_eo(x0, FuncAxStruct, ...
    FuncXStruct, Options)
% LBFGS - perform unconstrained minimization with the L-BFGS algorithm.


% Copyright 2010 Eli Osherovich.


%% Set parameters (default value)
[   x,...               % Initial x: obtained from x0 by converting
    ...                 % to column vector and changing type to
    ...                 % complex/real if necessary
    maxIter,...         % Maximal number of iterations (200)
    complexVarsFlag,... % Indicator whether the variables X are complex (false)
    nPrev,...           % Number of previous steps/gradients to remember (100)
    useMex,...          % Shall we use MEX files? (true)
    tolX, ...           % Step size tolerance (1e-8)
    tolFun, ...         % Function value tolerance (1e-8)
    tolGrad,...         % Gradient norm tolerance (1e-8)
    display...          % Progress report (false)
    ] = lbfgsGetOptions_eo(x0, Options);

%% Variables initialization and memory allocation.
% Preallocate space for previous steps and grads.
prevGrads = cell(1, nPrev);
prevSteps = cell(1, nPrev);

% Generate empty cell array of proper size (used by some functions).
empty = cell(size(FuncAxStruct));

% Calculate Ax for the initial point.
Ax = calculate_linop('forward', FuncAxStruct, x);

% Calculate initial function value and gradient.
[funcVal, grad] = calc_EDx(x, Ax, FuncAxStruct, FuncXStruct, empty, [], false, complexVarsFlag);
gradNorm = norm(grad);
funcCount = 1;


% Print initial state (if requested).
if display
	fprintf('%10s %10s %15s %15s %15s\n','Iteration','FunEvals','Step Length','Function Val','Opt Cond');
    fprintf('%10d %10d %15.5e %15.5e %15.5e\n', 0, ...
        funcCount, 0, funcVal, gradNorm);
end

% Test termination criteria.
[done, exitFlag, exitMsg] = testTermCriteria(...
    0, x, grad, funcVal, ...
    Inf, gradNorm, Inf,...
    [], maxIter, tolX, tolGrad, tolFun);

% Set output structure fileds.
Output.initialFval = funcVal;
Output.initialGradNorm = gradNorm;

% Initial step and search direction scale factor.
t0 =  1;
H0 = min(1, 1/norm(grad,1));

iter = 0;
validIdx = 0;
wrapAround = 0;
func =  @(x) minFunc_wrapper(x, FuncAxStruct, FuncXStruct);

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
    % Ad = calculate_linop('forward', FuncAxStruct, d);
    
    % Save current point data.
    gradOld = grad;
    funcValOld = funcVal;
    
    % Peform line search.
    [t, x, funcVal, grad, LSoutput] =  wolfeLS_eo(t0, x, funcVal, grad, ...
        func, d, Options);
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
    
    % Test termination criteria.
    [done, exitFlag, exitMsg] = testTermCriteria(...
        iter, x, grad, funcVal, ...
        stepNorm, gradNorm, abs(funcVal - funcValOld),...
        LSoutput, maxIter, tolX, tolGrad, tolFun);
   
end
% Display exit message (if requested).
if display
    fprintf('%s\n', exitMsg);
end

% Set output struct fields.
Output.finalFval = funcVal;
Output.firstOrderOpt = gradNorm;
Output.nIterations = iter;
Output.funcCount = funcCount;
Output.exitMsg = exitMsg;
Output.exitFlag = exitFlag;
Output.lsExitMsg = LSoutput.exitMsg;
Output.lsExitFlag = LSoutput.exitFlag; 

function [done, exitFlag, exitMsg] = testTermCriteria(...
    iter, x, grad, funcVal, ...
    stepNorm, gradNorm, funcValDiff,...
    LSoutput, maxIter, tolX, tolGrad, tolFun)

done = false;
exitFlag = 0;
exitMsg = '';

% Check if x has valid values.
if ~isValid(x)
    done = true;
    exitMsg = 'Invalid X values';
    exitFlag = -1;
    return;
end

% Check if grad has valid values.
if ~isValid(grad)
    done = true;
    exitMsg = 'Invalid Grad values';
    exitFlag = -2;
    return;
end

% Check if fVal is valid.
if ~isValid(funcVal)
    done = true;
    exitMsg = 'Invalid F value';
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
        exitMsg = sprintf('Line-search cound not find a suitable step-length');
        return;
    end
end
    

% If gradint's norm is small (below tolGrad) we probably found a local
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