function [maxIter, complexVarsFlag, display, saveDir] = ...
        hioGetOptions_eo(x0, options)


% shall we assume complex variables?
% by default use the "complexity" of x0.
complexVarsFlag = getOpt_eo(options, 'complexVarsFlag', ~isreal(x0));


% maximal number of iterations
maxIter = getOpt_eo(options, 'maxIter', 200);

% display type (progress report)
display = getOpt_eo(options, 'Display', true);

% allow also on/off settings for Display
if strcmpi(display, 'off')
    display = false;
elseif strcmpi(display, 'on')
    display = true;
end

% directory to save current x (empty by default = not to save)
saveDir = getOpt_eo(options, 'saveDir', '');
% create the directory if it does not exist
if ~isempty(saveDir) && ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
