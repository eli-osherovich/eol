function [t, x, f, g, Ax, Output] = wolfeLS_eo(t, x, f, g, dir, Ax, Adir, FuncAxStruct, funcX, Options)
% Line-search routine with strong Wolfe's conditions.
    
    
    
% Copyright 2010 Eli Osherovich.
    
    
% Save initial data (at t=0).
t0 = 0;
x0 = x;
f0 = f;
g0 = g;
d = real(g0'*dir); % directional derivative
d0 = d; % keep the directional derivative at the starting point
Ax0 = Ax;    
    
% Get parameters.
[   maxIter,...         % Maximal number of iterations (20)
    fdf, ...            % Function value decrease constant (1e-4)
    gdf, ...            % Gradient norm decrease constant (0.9)
    tolX, ...           % Interval of uncertainty width (1e-8)
    tolDD, ...          % Directional derivative tolerance (1e-10)
    minFVal, ...        % Minimal function value (-Inf)
    minStp, ...         % Minimal step (1e-20)
    maxStp, ...         % Maximal step (1e+20 or based on minFVal and
    xtrapMin, ...       % Step-length extrapolation max. factor (1.1)
    xtrapMax, ...       % Step-length extrapolation min. factor (20)
    complexVarsFlag...  % Indicator whether the variables X are complex (false)
] = wolfeLSGetOptions_eo(x0, Options);





% We shall keep all function/gradient values obtained during this
% subroutine: each row of data contains a triplet (t, f, d).
data = zeros(maxIter+1, 3);
data(1, :) = [t0, f0, d0];

% Update maxStp if minimal function value was given.
maxStp = min(maxStp, -(f0 - minFVal)/(d0*fdf));


% Check inputs validity.
if t <=0 || fdf < 0 || gdf < 0 || tolX < 0 || tolDD < 0 || ...
       minStp < 0 || maxStp < minStp || maxIter < 0
    if nargout < 5
        error('EOL:wolfeLS:InvalidInputs', 'Invalid input parameters');
    else
        Output.exitFlag = -2;
        Output.exitMsg = 'Invalid input parameters';
        Output.funcCount = 0;
        Output.nIterations = 0;
        return;
    end
end

% Set initial X_lo data.
idx_lo = 1; t_lo = t0; f_lo = f0;

% Bracket is not found hence i_hi is empty.
idx_hi = [];

iter = 0;
done = false;

%% Run iterations
while ~done
    
    % Force the step to be within global bounds.
    t = max(t, minStp);
    t = min(t, maxStp);
    
    
    % Test termination criteria.
    [done, exitFlag, exitMsg] = testTermCriteriaLS(...
        data, iter, idx_lo, idx_hi, f, d, ...
        maxIter, tolX, tolDD, minFVal, maxStp, minStp);
    
    % Update iteration counter and index of the current slot in DATA array. 
    iter = iter + 1;
    idx = iter + 1;
    
    % Save current *_lo data.
    idx_lo_old = idx_lo;
    t_lo_old = t_lo;
    
    % Evaluate the function and gradient at the new point.
    x = x0 + t*dir;
    %     for k = 1:length(FuncAxStruct)
    %             Ax{k} = Ax0{k} + t * Adir{k};
    %     end
    Ax = applyMapping(FuncAxStruct, x);
    [f, g] = calcObjFunc(x, Ax, FuncAxStruct, funcX, complexVarsFlag);
    
    d = real(g'*dir);
    
    % Save current data.
    data(idx, :) = [t, f, d];
    
    if f > f0 + fdf*t*d0 || f >= f_lo
        % Sufficient decrease condition is not statisfied
        % or current function value is higher than f_lo.
        % Bracket is found.
        idx_hi = idx;
        
    elseif abs(d) <= abs(gdf*d0)
        % Current t satisfies sufficient decrease condition
        % if it also satisfies curvature condition. Hence, the strong
        % Wolfe conditions are satisfied. 
        % Line-search succeeded.
        break;
    else
        % Current t satisfies the sufficient decrease condition and
        % corresponding function value is below f_lo. Hence, *_lo
        % variables must be updated.
        idx_lo = idx; t_lo = t; f_lo = f;
        
        % Check if we have a bracket here.
        if d*(t_lo_old - t) < 0
            idx_hi = idx_lo_old;
        end
        
    end
    
    % Set the maximal and minimal step length for the current iteration.
    if ~isempty(idx_hi)
        % If bracket is known, use its endpoints, ...
        tMinIter = min(data(idx_lo,1), data(idx_hi, 1));
        tMaxIter = max(data(idx_lo,1), data(idx_hi, 1));
    else
        % ... otherwise, use some extrapolated values.
        tMinIter = t + xtrapMin * (t - t_lo_old);
        tMaxIter = t + xtrapMax * (t - t_lo_old);
    end
    
    t = step_eo(data, idx, idx_lo, idx_hi, tMinIter, tMaxIter);
end


Output.finalFval = f;
Output.firstOrderOpt = abs(d);
Output.nIterations = iter;
Output.funcCount = iter;
Output.exitFlag = exitFlag;
Output.exitMsg = exitMsg;

function [done, exitFlag, exitMsg] = testTermCriteriaLS(...
    data, iter, i_lo, i_hi, f, d, ...
    maxIter, tolX, tolDD, minFVal, maxStp, minStp)
    
    done = false;
    exitFlag = 0;
    exitMsg = '';
    
    t = data(iter+1, 1);
    
    if iter == 0 && d >= 0
        done = true;
        exitFlag = -1;
        exitMsg = 'Directional derivative is not a descent direction';
        return;
    end
    
    if iter >= maxIter
        done = true;
        exitFlag = 1;
        exitMsg = sprintf(['Number of iterations (%d) exceeded maximum ' ...
                           'allowed.'], maxIter);
        return;
    end
    
    if ~isempty(i_hi)
        t_lo = data(i_lo, 1);
        d_lo = data(i_lo, 3);
        t_hi = data(i_hi, 1);
        
        % Check that the bracket is valid, i.e., contains minimum.
        if d_lo * ( t_hi - t_lo) > 0 
            done = true;
            exitFlag = -4;
            exitMsg = 'Invalid bracket';
            return;
        end
        
        % Check that the bracket is not too short.
        if abs(t_lo - t_hi) <= tolX
            done = true;
            exitFlag = 2;
            exitMsg = sprintf(['Interval of uncertainty is shorter than ' ...
                               'tolX (%g) tolerance.'], tolX);
            return;
        end
    end
    
    if f <= minFVal
        done = true;
        exitFlag = 3;
        exitMsg = sprintf('Function value is below minFVal (%g).', ...
                          minFVal);
        return;
    end
    
    if t == maxStp
        done = true;
        exitFlag = 4;
        exitMsg = sprintf('t is equal maxStp (%g)', maxStp);
        return;
    end
    
    if t == minStp
        done = true;
        exitFlag = 5;
        exitMsg = sprintf('t is equal minStp (%g)', minStp);
        return;
    end
    
    if abs(d) <= tolDD
        done = true;
        exitFlag = 6;
        exitMsg = sprintf(['Directional derivative is below tolDD (%g) ' ...
                           'tolerance'], tolDD);
        return;
    end
    
