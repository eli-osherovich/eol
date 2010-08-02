function [x, fval, exitflag, output] = lbfgsO_eo(x0, func_Ax_Struct, ...
                                                func_x_Struct, options)
% LBFGS - perform unconstrained minimization with the L-BFGS algorithm.


% Copyright 2010 Eli Osherovich.


%% set parameters (default value)
[   x,...               % Initial x: obtained from x0 by converting
    ...                 % to column vector and changing type to
    ...                 % complex/real if necessary
    maxIter,...         % Maximal number of iterations (200)
    ComplexVarsFlag,... % Indicator whether the variables X are complex (false)
    nPrev,...           % Number of previous steps/gradients to remember (100)
    useMex,...          % Shall we use MEX files? (true)
    tolX, ...           % Step size tolerance (1e-8)
    tolFun, ...         % Function value tolerance (1e-8)
    tolGrad,...         % Gradient norm tolerance (1e-8)
    display...          % Progress report (false)
    ] = lbfgsGetOptions_eo(x0, options);

%% Variables initialization and memory allocation
% previous steps and grads
PrevGrads = cell(1, nPrev);
PrevSteps = cell(1, nPrev);

% generate empty cell array of proper size (used by some functions)
empty = cell(size(func_Ax_Struct));

% calculate Ax for the initial point
Ax = calculate_linop('forward', func_Ax_Struct, x);

% calculate initial function value and gradient
[fval, grad] = calc_EDx(x, Ax, func_Ax_Struct, func_x_Struct, empty, [], false, ComplexVarsFlag);
grad_norm = norm(grad);
funcCount = 1;


% print initial state (if requested)
if display
	fprintf('%10s %10s %15s %15s %15s\n','Iteration','FunEvals','Step Length','Function Val','Opt Cond');
    fprintf('%10d %10d %15.5e %15.5e %15.5e\n', 0, ...
        funcCount, 0, fval, grad_norm);
end

% test termination criteria
[done, exitflag, exitmsg] = testTermCriteria(...
    0, Inf, grad_norm, Inf, maxIter, tolX, tolGrad, tolFun);

% set output structure fileds
output.initialFval = fval;
output.initialGradNorm = grad_norm;

t0 =  1;
H0 = min(1,1/sum(abs(grad)));

iter = 0;
validIdx = 0;
wrapAround = 0;
func =  @(x) minFunc_wrapper(x, func_Ax_Struct, func_x_Struct);

%% Run iterations
while ~done
    % update iteration counter
    iter = iter + 1;
    
    % calculate new search direction
    if useMex
        d = lbfgsDirPersistC_eo(validIdx, wrapAround, -grad, ...
            PrevSteps, PrevGrads, H0);
    else
        d = lbfgsDirPersist_eo(validIdx, wrapAround, -grad, ...
            PrevSteps, PrevGrads, H0);
    end
    
    % calculate d's linear transform (used by line search)
    Ad = calculate_linop('forward', func_Ax_Struct, d);
    
    % save current point data
    %x_old = x;
    grad_old = grad;
    fval_old = fval;
    
    % peform line search
   
    [t,fval,grad,LSfunEvals] = WolfeLineSearch_mf(x,t0,d,fval,grad,grad'*d,1e-4,0.9,4,25,tolX,0,0,1,func);
    x = x+t*d;
    funcCount = funcCount + LSfunEvals;
    
     
    % go to the new point
    grad_norm = norm(grad);
    step = t*d;
    step_norm = norm(step);
    
    
    % update LBFGS memory
    [PrevSteps,PrevGrads,H0, validIdx, wrapAround] = ...
        lbfgsUpdate_eo(validIdx, wrapAround, step, grad-grad_old,...
        PrevSteps, PrevGrads,H0);

    % print progress (if requested)
    if display
        fprintf('%10d %10d %15.5e %15.5e %15.5e\n', iter, ...
            funcCount, t, fval, grad_norm);
    end
    
    % test termination criteria
    [done, exitflag, exitmsg] = testTermCriteria(...
        iter, step_norm, grad_norm, abs(fval - fval_old),...
        maxIter, tolX, tolGrad, tolFun);
end
% display exit message (if requested)
if display
    fprintf('%s\n', exitmsg);
end

% set output struct fields
output.finalFval = fval;
output.firstorderopt = grad_norm;
output.iterations = iter;
output.funcCount = funcCount;
output.exitmsg = exitmsg;
output.exitflag = exitflag;

function [done, exitflag, msg] = testTermCriteria(...
    iter, stepNorm, gradNorm, fvalDiff, maxIter, tolX, tolGrad, tolFun)

done = false;
exitflag = 0;
msg = '';

if (gradNorm < tolGrad)
    done = true;
    exitflag = 1;
    msg = sprintf(['Gradient magnitude is smaller that tolGrad (%g) ' ...
                   'tolerance.'], tolGrad);
    return;
end

if (stepNorm < tolX)
    done = true;
    exitflag = 2;
    msg = sprintf('Step length is shorter than tolX (%g) tolerance.', ...
                  tolX);
    return;
end

if (fvalDiff < tolFun)
    done = true;
    exitflag = 3;
    msg = sprintf(['Change in the obejctive function value is smaller ' ...
                   'than tolFun (%g) tolerance.'], tolFun);
    return;
end

if (iter > maxIter)
    done = true;
    exitflag = 0;
    msg = sprintf(['Number of iterations (%d) exceeded maximum ' ...
                   'allowed.'], maxIter);
    return;
end
