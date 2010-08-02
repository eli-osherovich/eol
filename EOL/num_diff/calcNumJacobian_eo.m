function Jt = calcNumJacobian_eo(func, x, mode)
% CALCNUMJACOBIAN - calculate Jacobian numerically.
%
% Jt =  CALCNUMJACOBIAN(FUNC, X) calculate Jacobian matrix numerically by the
% central difference approximation. 
%
% The output JT is the conjugate transpose of the Jacobian (gradient)
% matrix of FUNC computed at point X.


% Copyright 2007-2010 Eli Osherovich.


% Use central difference by default.
if nargin < 3
    mode = 'central';
end

% Number of variables.
N = numel(x);

% Calculate func(x) for current x.
% The result is used to determine the number of equations (M) and in
% forward/backward mode.
fx = func(x);

% Get number of equations.
M = numel(fx);
% Preallocate space for Jt.
Jt = zeros(N, M);
    
    
% Make mode case insensitive.
mode = lower(mode);

% Set step size.
h = findOptH_eo(x, mode);

% Perform actual calculations.
switch mode
  case {'forward', 'backward'}
    % Calculate partial derivatives df/dx_i (forward/backward difference).
    for i = 1:N
        xTmp = x;
        if mode(1) == 'f'
            % Forward difference.
            xTmp(i) = x(i) + h(i);
            Jt(i,:) = (func(xTmp) - fx)/h(i);
        else
            % Backward difference.
            xTmp(i) = x(i) - h(i);
            Jt(i,:) = (fx - func(xTmp))/h(i);
        end
    end
    
  case 'central'
    % Calculate partial derivatives df/dx_i (central difference).
    for i = 1:N
        xF = x;
        xB = x;
        xF(i) = x(i) + h(i);
        xB(i) = x(i) - h(i);
        Jt(i,:) = (func(xF) - func(xB))/2/h(i);
    end
    
  case 'precise'
    % Calculate partial derivatives df/dx_i (using central difference
    % with Richardson's extrapolation).
    
    % Maximal number of extrapolation steps.
    nSteps = 15;
    
    % Preallocate space for all Jacobians.
    Jall = zeros(N, M, nSteps);
    
    % start with the maximal step
    h = h * 2^nSteps;
    
    % Calculate Jacobians for all step sizes.
    for s = 1:nSteps
        % Calculate partial derivatives df/dx_i (central difference).
        for i = 1:N
            xF = x;
            xB = x;
            xF(i) = x(i) + h(i);
            xB(i) = x(i) - h(i);
            Jall(i,:,s) = (func(xF) - func(xB))/2/h(i);
        end
        h = h/2;
    end
    
    % Jall now contains all possible Jacobians (with all step sizes).
    % Perform Richardson's extrapolation.
    Jt = extRichardson(Jall);
    
  otherwise
     error('EOL:CALCNUMJACOBIAN:WrongInput', ['Unknown mode: %s it must ' ...
                         'be either forward, backward, or central'], mode);
end


function Jt = extRichardson(Jall)

nSteps = size(Jall, 3);
Jsize = [size(Jall, 1), size(Jall,2)];
corrections = cell(nSteps - 1);
Jtableau = cell(nSteps);



minErr = Inf(Jsize);
Jtableau{1,1} = Jall(:,:,1);
Jt = Jall(:,:,1);

for m = 2:nSteps
    Jtableau{m,1} = Jall(:,:,m);
    good_idx = true(Jsize);

    p_i = 2;
    for n = 2:m
        frac = 2^p_i - 1;
        
        Jtableau{m,n} = ((2^p_i)*Jtableau{m, n-1} - Jtableau{m-1, n-1})/frac;
        corrections{m-1, n-1} = Jtableau{m,n} - Jtableau{m, n-1};
    
        if n > 2
            good_idx(abs(corrections{m-1,n-1}) >= abs(corrections{m-1,n-2})) = false;
        end
        errDecreaseIdx = minErr > abs(corrections{m-1,n-1});
        Jt(good_idx & errDecreaseIdx) = Jtableau{m,n}(good_idx & errDecreaseIdx);
        minErr(errDecreaseIdx) = abs(corrections{m-1,n-1}(errDecreaseIdx));
        
        p_i = p_i + 2;
    end
    
end
