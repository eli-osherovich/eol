function [...
    x0, maxIter, complexVarsFlag, nPrev, useMex,...
    tolX, tolFun, tolGrad, display] = lbfgsGetOptions_eo(x0, options)

% convert x to a column vector and verify validity

if isnumeric(x0)
    if ~isa(x0, 'double')
        warning('EOL:LBFGS:NonDoubleInput',...
            'LBFGS supports only double variables. Converting to double.');
        x0 = double(x0(:));
    else
        x0 = x0(:);
    end
else
    error('EOL:LBFGS:NonNumericInput',...
        'X0 must be of a numeric type.')
end

% shall we assume complex variables?
% by default use the complexity of x0
complexVarsFlag = getOpt_eo(options, 'complexVarsFlag', ~isreal(x0));

% check whether the type of x matches complexVarsFlag
if complexVarsFlag && isreal(x0)
    % complexVarsFlag is set explicitly but x0 is real
    warning('EOL:LBFGS:RealXwithcomplexVarsFlag',...
        'complexVarsFlag is set while X0 is real. Converting to complex.');
    x0 = complex(x0, zeros(size(x0)));
end

% maximal number of iterations
maxIter = getOpt_eo(options, 'maxIter', 200);

% number of previous steps/gradients to use
nPrev = getOpt_eo(options, 'nPrev', 100);

% shall MEX files be used
useMex = getOpt_eo(options, 'useMex', true);

% step norm tolerance
tolX = getOpt_eo(options, 'tolX', 1e-8);

% function change tolerance
tolFun = getOpt_eo(options, 'tolFun', 1e-8);

% gradient norm tolerance
tolGrad = getOpt_eo(options, 'tolGrad', 1e-8);

% display type (progress report)
display = getOpt_eo(options, 'Display', true);

% allow also on/off settings for Display
if strcmpi(display, 'off')
    display = false;
elseif strcmpi(display, 'on')
    display = true;
end
