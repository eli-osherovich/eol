function [x, err, Output] = HIO_eo(x0, supportX, FMod, supportF, phaseLo, phaseHi,...
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
    
   % Check if the third output argument (Output) was requested.
   if nargout > 2
       outputRequested = true;
       Output.allFourierErrs = NaN(maxIter+1, 1);
       Output.allObjectErrs = NaN(maxIter+1, 1);
   else
       outputRequested = false;
   end

    
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
    % Print initial state (if requested).
    if display
        fprintf('%10s %15s %15s %15s %15s\n','Iteration', 'Step Length','Total Err','Err Fourier', 'Err Signal');
    end
    
    for i = 1:maxIter
        
        % Generate singnal that satisfies Fourier domain constraints.
        [y, errF] = forceFourierMagnitudeAndPhase(x, xSize, FMod, supportF, ...
            complexVarsFlag, phaseLo, phaseHi);
        
        
        % Support violations. We use abs(y) since y can be complex.
        violations = ~supportX(:) & abs(y);
        
        % Non-negativity violations (for real signals).
        if ~complexVarsFlag
            assert(isreal(y));
            violations = violations | (supportX(:) & y<0);
        end
        
        % Signal domain error.
        errS = norm(y(violations))^2;
        
        % Total error.
        err = errF + errS;
    
        % If Output structure was requested update the list of all errors.
        if outputRequested
            Output.allFourierErrs(i + 1) = errF;
            Output.allObjectErrs(i + 1) = errS;
        end
        
        % Print progress (if requested).
        if display
            fprintf('%10d %15.5e %15.5e %15.5e %15.5e\n', i, ...
                0, err, errF, errS);
        end
        
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
