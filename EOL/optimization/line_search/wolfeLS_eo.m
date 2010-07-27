function [t, x, f, g, Output] = wolfeLS_eo(t, x, f, g, ...
    func, dir, Options)
% Line-search routine with strong Wolfe's conditions.



% Copyright 2010 Eli Osherovich.


% Save initial data (at t=0).
t0 = 0;
x0 = x;
f0 = f;
g0 = g;
d = real(g0'*dir); % directional derivative
d0 = d;


% Get parameters.
[   maxIter,...         % Maximal number of iterations (20)
    fdf, ...            % Function value decrease constant (1e-4)
    gdf, ...            % Gradient norm decrease constant (0.9)
    tolX, ...           % Interval of uncertainty width (1e-8)
    tolDD, ...          % Directional derivative tolerance (1e-10)
    minFVal, ...        % Minimal function value (-Inf)
    minStp, ...         % Minimal step (1e-20)
    maxStp, ...         % Maximal step (1e+20 or based on minFval and
    xtrapMin, ...       % Step-length extrapolation max. factor (10)
    xtrapMax, ...       % Step-length extrapolation min. factor (1.1)
    complexVarsFlag...  % Indicator whether the variables X are complex (false)
    ] = wolfeLSGetOptions_eo(x0, Options);
 



% dummy empty cell array (used by calc_EDx)    
% empty = cell(size(func_Ax_Struct));

% pre-allocate space for Ax
% Ax = cell(size(Ax0));



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
        return;
    end
end

% Set initial X_lo data.
i_lo = 1; t_lo = t0; f_lo = f0;

% Bracket is not found hence i_hi is empty.
i_hi = [];

iter = 0;
done = false;

%% Run iterations
while ~done
    % Test termination criteria.
    [done, exitFlag, exitMsg] = testTermCriteriaLS(...
        iter, data, i_lo, i_hi, f, d, ...
        maxIter, tolX, tolDD, minFVal, maxStp, minStp);
    
    
    iter = iter + 1;
    
    % Save current X_lo data.
    i_lo_old = i_lo;
    t_lo_old = t_lo;
    
    % Force the step to be within global bounds.
    t = max(t, minStp);
    t = min(t, maxStp);
        
    
        
    % evaluate function and gradient at the new point
    x = x0 + t*dir;
    %     for k = 1:length(func_Ax_Struct)
    %         Ax{k} = Ax0{k} + t * Adir{k};
    %     end
    %     [f, g] = calc_EDx(x, Ax, func_Ax_Struct, func_x_Struct, empty, [], false, ComplexVarsFlag);
    
    [f, g] = func(x);
    d = real(g'*dir);
    
    % Save current data.
    data(iter+1, :) = [t, f, d];
    
    if f > f0 + fdf*t*d0 || f > f_lo
        % sufficient decrease condition is not statisfied
        % or current function value is higher than f_lo
        % bracket is found
        i_hi = iter+1;
        
    elseif abs(d) <= abs(gdf*d0)
        % current t satisfies sufficient decrease condition
        % if it also satisfies curvature condition it satisfies strong
        % Wolfe criteria. line-search succeeded.
       break;
    else
        % current t satisfies sufficient decrease conditions and current
        % function value is below f_lo. hence, f_lo must be updated
        i_lo = iter+1; t_lo = t; f_lo = f;
        
        if d*(t_lo_old - t) < 0
            i_hi = i_lo_old;
        end
        
    end
    
    % Set the maximal and minimal step length for the current iteration.
    if ~isempty(i_hi)
        % If bracket is known, use its endpoints, ...
        tMinIter = min(data(i_lo,1), data(i_hi, 1));
        tMaxIter = max(data(i_lo,1), data(i_hi, 1));
    else
        % ... otherwise, use some extrapolated values.
        tMinIter = t + xtrapMin * (t - t_lo_old);
        tMaxIter = t + xtrapMax * (t - t_lo_old);
    end
    
    t = step_eo(iter+1, data, i_lo_old, i_hi, tMinIter, tMaxIter);
end


Output.finalFval = f;
Output.firstOrderOpt = abs(d);
Output.nIterations = iter;
Output.funcCount = iter;
Output.exitFlag = exitFlag;
Output.exitMsg = exitMsg;

function [done, exitFlag, exitMsg] = testTermCriteriaLS(...
    iter, data, i_lo, i_hi, f, d, ...
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
        exitMsg = sprintf('Interval of uncertainty is shorter than tolX (%g) tolerance.', ...
            tolX);
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
    exitMsg = sprintf('Directional derivative is below tolDD (%g) tolerance', ...
        tolDD);
    return;
end
