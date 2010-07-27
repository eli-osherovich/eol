function [maxIter, fdc, gdf, tolX, tolDD, minFVal, minStp, maxStp,...
        xtrapMin, xtrapMax, complexVarsFlag] = wolfeLSGetOptions_eo(x0, Options)
% WOLFELSGETOPTIONS - get options for wolfeLS routine.


% Copyright 2010 Eli Osherovich.
    
    
% Maximal number of iterations.
maxIter = getOpt_eo(Options, 'LSmaxIter', 20);

% (Sufficient) Function decrease constant.
fdc = getOpt_eo(Options, 'LSFvalDC', 1e-4);

% (Sufficient) Gradient norm decrease factor.
gdf = getOpt_eo(Options, 'LSGradDF', 0.9);

% Interval of uncertainty tolerance.
tolX = min(...
    getOpt_eo(Options, 'LSTolX', Inf),...
    getOpt_eo(Options, 'tolX', 1e-8));

% Directional derivative tolerance.
tolDD = getOpt_eo(Options, 'LSTolDD', 1e-10);

% Minimum function value. Should be set if known. 
% For example, in the least squares problems minFVal = 0.
minFVal = getOpt_eo(Options, 'minFVal', -Inf);

% minimal step
minStp = getOpt_eo(Options, 'LSMinStep', 1e-20);

% maximal step
maxStp = getOpt_eo(Options, 'LSMaxStep', 1e+20);

% min. extrapolation factor
xtrapMin = getOpt_eo(Options, 'LSXtrapMin', 1.1);

% max. extrapolation factor
xtrapMax = getOpt_eo(Options, 'LSXtrapMax', 10);

% shall we assume complex variables?
% by default use the complexity of x0
complexVarsFlag = getOpt_eo(Options, 'complexVarsFlag', ~isreal(x0));
