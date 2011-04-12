function [x, err] = HIO_eo(x0, supportX, FMod, supportF, phaseLo, phaseHi,...
        Options)
    %HIO - hybrid input-output algorithm with padding

    
    
    % Copyright 2008-2011 Eli Osherovich.
    
    
    
           
    % Notation:
    % x - (current) signal estimate.
    % y - transformation (projection) of x that satisfies Fourier domain
    % constraints.
    
    
    % Check inputs.
    validateattributes(x0, {'numeric'},{'nonnan', 'finite'});
    validateattributes(supportX, {'logical'},{});
    validateattributes(FMod, {'numeric'}, {'nonnegative','nonnan', 'finite'});
    validateattributes(supportF, {'logical'},{});
    

    
    %% Use default options if Options is not provided or empty.
    if nargin < 7
        Options = struct();
    end
    
    %% Set parameters (default value)
    [   maxIter,...         % Maximal number of iterations (1000)
        complexVarsFlag,... % Indicator whether the variables X are complex (false)
        display,...         % Progress report (true)
        saveDir...          % Directory to save current x ('' = do not save)
    ] = hioGetOptions_eo(x0, Options);
    
    
    
    % Damping factor, several publications claim that it must be around
    % .75
    beta = 0.75;

    xSize = size(x0);
    
    
    % Run HIO iterations.
    x = x0;
    % Save current x (if requested).
    if ~isempty(saveDir)
        save(fullfile(saveDir, int2str(0)), 'x');
    end
    
    for i = 1:maxIter
        
        % Generate singnal that satisfies Fourier domain constraints.
        [y, errF] = forceFourierMagnitudeAndPhase(x, xSize, FMod, supportF, ...
            complexVarsFlag, phaseLo, phaseHi);
        
        
        % Support violations. We use abs(y) since y can be complex.
        violations = ~supportX(:) & abs(y);
        
        % Non-negativity violations (for real signals).
        if ~complexVarsFlag
            violations = violations | (supportX(:) & y<0);
        end
        
        % Signal domain error.
        errS = norm(y(violations));
        
        % Total error.
        err = errF + errS;
        
        % HIO step.
        y(violations) = x(violations) - beta*y(violations);
        x = y;
        
        
        % Save current x (if requested).
        if ~isempty(saveDir)
            save(fullfile(saveDir, int2str(i)), 'x');
        end
    end
    
    % Reshape x to the original size.
    x = reshape(x, size(x0));
