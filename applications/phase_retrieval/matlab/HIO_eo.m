function [x, err] = HIO_eo(x0, supportX, FMod, supportF, nIter,...
        complexVarsFlag, phaseLo, phaseHi)
    %HIO - hybrid input-output algorithm with padding

    
    
    % Copyright 2008-2011 Eli Osherovich.
    
    
    
           
    % Notation:
    % x - (current) signal estimate.
    % y - transformation (projection) of x that satisfies Fourier domain
    % constraints.
   
    
    % Check inputs.
    validateattributes(x0, {'numeric'},{});
    validateattributes(supportX, {'logical'},{});
    validateattributes(FMod, {'numeric'}, {'nonnegative'});
    validateattributes(supportF, {'logical'},{});
    validateattributes(complexVarsFlag, {'logical'}, {'scalar'});
    if nargin ~= 6 && nargin ~=8
        error('You must provide either 6 or 8 inputs');
    end
    
    if nargin ~= 8
        phaseLo = [];
        phaseHi = [];
    end
    
    
    % Damping factor, several publications claim that it must be around
    % .75
    beta = 0.75;

    xSize = size(x0);
    
    
    % Run HIO iterations.
    x = x0;
    for i = 1:nIter
        
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
    end
    
    % Reshape x to the original size.
    x = reshape(x, size(x0));
