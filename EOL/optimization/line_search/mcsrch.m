function [t, x, f, g, output] =  mcsrch(t, x0, f0, g0, ...
        func, d, options)
%MCSRCH - More & Thuente line search routine

% get parameters
[   maxIter,...         % Maximal number of iterations (20)
    fdf, ...            % Function value decrease constant (1e-4)
    gdf, ...            % Gradient norm decrease constant (0.9)
    tolX, ...           % Interval of uncertainty width (1e-8)
    minFval, ...        % Minimal function value (-Inf)
    minStp, ...         % Minimal step (1e-20)
    maxStp, ...         % Maximal step (1e+20 or based on minFval and
    xtrapF, ...         % Step-length extrapolation factor (10)
    complexVarsFlag...  % Indicator whether the variables X are complex (false)
    ] = wolfeLSGetOptions_eo(x0, options);


if ~libisloaded('mcsrch')
    loadlibrary('line_search/mcsrch.so', 'mcsrch.h');
end

% prepare inputs
n = libpointer('int32Ptr', numel(x0));
x = libpointer('doublePtr', x0);
f = libpointer('doublePtr', f0);
g = libpointer('doublePtr', g0);
s = libpointer('doublePtr', d);
stp = libpointer('doublePtr', t);
ftol = libpointer('doublePtr', fdf);
xtol = libpointer('doublePtr', tolX);
maxfev = libpointer('int32Ptr', maxIter);
info = libpointer('int32Ptr', 0);
nfev = libpointer('int32Ptr', 0);
gtol = libpointer('doublePtr', gdf);
stpmin = libpointer('doublePtr', minStp);
stpmax = libpointer('doublePtr', maxStp);
wa = libpointer('doublePtr', zeros(numel(x0), 1));

while true
    calllib('mcsrch', 'mcsrch_', n, x, f, g, s, stp, ftol, xtol, ...
        maxfev, info, nfev, gtol, stpmin, stpmax, wa);
    
    switch info.value
        case -1
            % a return is made to compute the function and gradient
            [f, g] = func(x.value);
            f = libpointer('doublePtr', f);
            g = libpointer('doublePtr', g);
        case 0
            % improper input parameters
            error('EOL:mcsrch:InvalidInputs', 'Invalid input paramters');
        case 1
            % sufficient decrease condition and the directional derivative
            % condition hold
            break
        case {2, 3, 4, 5, 6}
            % 2: relative width of the interval of uncertainty is at most xtol
            % 3: number of calls to func has reached maxev
            % 4: the step is at the lower bound minStp
            % 5: the step is at the upper bound maxStp
            % 6: rounding errors prevent further progress. there may not be a
            %    step which satisfies the sufficient decrease and curvature
            %    conditions. tolerance may be too small.
            break
        case 7
            % d is not a descent direction
            error('EOL:mcsrch:WrongDirection', 'Not a descent direction');
        otherwise
            % should never get here
            error('EOL:mcsrch:UnknownError', 'Unknown error ocurred');
    end
end

            

t = stp.value;
x = x.value;
f = f.value;
g = g.value;


output.finalFval = f;
output.firstorderopt = real(g'*d);
output.iterations = nfev.value;
output.funcCount = nfev.value;
output.exitmsg = '';
output.exitflag = info.value;
