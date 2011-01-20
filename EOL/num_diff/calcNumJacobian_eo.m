function J = calcNumJacobian_eo(x, func, mode)
% CALCNUMJACOBIAN - calculate Jacobian numerically.
%
% Jt =  CALCNUMJACOBIAN(X, FUNC) calculate Jacobian matrix numerically by the
% central difference approximation. 
%
%
% Assumptions:
% -------------
% x  must be a column vector


% Copyright 2007-2011 Eli Osherovich.




% Make sure x is a numeric vector
validateattributes(x, {'numeric'}, {'vector'});

% Use central difference by default.
if nargin < 3
    mode = 'central';
end






    
    
% Make mode case insensitive.
mode = lower(mode);

% Determine optimal/initial step size (done separately for the real and
% imaginary parts)
if isreal(x)
    hRe = findOptH_eo(x, mode);
    hIm = [];
else
    hRe = findOptH_eo(real(x), mode);
    hIm = 1i * findOptH_eo(imag(x), mode);
end

% Perform actual calculations.
% Note that hRe and hIm are assumed to be column vectors (hIm can be
% empty).
Jt = 0;
for h = [hRe hIm]
    switch mode
        case {'forward'}
            % Forward difference
            Jt = Jt + onesideDerivative(x, func, h);
        
        case {'backward'}
            % Backward difference
            Jt = Jt + onesideDerivative(x, func, -h);
        
        case 'central'
            Jt = Jt + centralDerivative(x, func, h);
            
        case 'precise'
            Jt = Jt + extrapDerivative(x, func, h);
            
            
        otherwise
            error('EOL:CALCNUMJACOBIAN:WrongInput', ['Unknown mode: %s it must ' ...
                'be either forward, backward, or central'], mode);
    end
end
J = Jt.';

function Jt = onesideDerivative(x, func, h)
% Number of variables.
N = numel(x);

% Calculate func(x) for current x.
fx = func(x);

% Get number of equations.
M = numel(fx);

% Preallocate space for Jt.
Jt = zeros(N, M);

for i = 1:N
    xTmp = x;
    
    % Change current x's entry.
    xTmp(i) = x(i) + h(i);
    
    % Calculate numerical derivative.
    Jt(i,:) = (func(xTmp) - fx)/h(i);
end

function Jt = centralDerivative(x, func, h)
% Number of variables.
N = numel(x);

% We need to get the number of equations (M). To avoid unnecessary function
% calculations we compute the first column of Jt separately from the rest.
xF = x; xF(1) = x(1) + h(1);
xB = x; xB(1) = x(1) - h(1);
J1 = (func(xF) - func(xB))/2/h(1);
M = numel(J1);

% Preallocate space for Jt and fill the first column.
Jt = zeros(N, M);
Jt(1,:) = J1;

% Get number of equations.
for i = 2:N
    xF = x;
    xB = x;
    xF(i) = x(i) + h(i);
    xB(i) = x(i) - h(i);
    % See note above (forward difference).
    Jt(i,:) = (func(xF) - func(xB))/2/h(i);
end

function Jt = extrapDerivative(x, func, h)
% Calculate derivative approximation using central difference and
% Richardson's extrapolation.
%
% Organization details:
% ---------------------
% The data is orginized in two tabeleaus: J (for the Jacobian values) and C
% (for Richardson's corrections). One can visualize these two as a single 
% tableau as shown below.
%
% [j11                                                ]
% [j21 (+) c11 (=) j22                                ]
% [j31 (+) c21 (=) j32 (+) c22 (=) j33                ]
% [j41 (+) c32 (=) j42 (+) c32 (=) j43 (+) c33 (=) j44]
%


% Maximal number of extrapolation steps.
nSteps = 15;

% Tolerance for the corrections quotient (along a column).
qTol = 2;


% Start with the maximal step.
h = h * 2^nSteps;

% Preallocate space for tableaus
Jtab= cell(nSteps);
Ctab = cell(nSteps-1);

% The first approximation is computed outside the look just to give us the
% size of the Jacobian.
Jtab{1,1} = centralDerivative(x, func, h);
h = h/2;

Jt = Jtab{1,1};
Jsize = size(Jt);

% Error estimate (separate for each entry of the Jacobian matrix).
minErr = inf(Jsize);




for i = 2:nSteps
    % Keep the index of "good" corrections:
    % 1) corrections must be monotonically decreasing along a row
    % 2) their quotient along a column must be 2^p (for the current p)
    cGoodIdx = true(Jsize);
        
    % Calculate the Jacobian corresponding to the current step.
    Jtab{i,1} = centralDerivative(x, func, h);
    
    % Calculate (the rest of) the current row.
    % In the central difference approximation the power series contains all
    % even powers starting from 2. Hence we we the current power to 2 and
    % then update it for every column
    p = 2;
    for j = 2:i
        Ctab{i-1, j-1} = (Jtab{i, j-1} - Jtab{i-1, j-1})/(2^p - 1);
        Jtab{i, j} = Jtab{i, j-1} +  Ctab{i-1, j-1};
        
        
        % Corrections that are not monotonically decreasing marked as 
        % "not good"
        if j > 2
            cGoodIdx(abs(Ctab{i-1, j-1}) > abs(Ctab{i-1, j-2})) = false;
        end
        % Corrections that does not obey the quotient (along the column) of
        % 2^p are marked as "not good"
        if i > j
            q = abs(Ctab{i-2, j-1}) ./ abs(Ctab{i-1, j-1});
            cGoodIdx(q > 2^p * qTol | q < 2^p / qTol) = false;
        end
        
        % Find corrections that provide new (better=smaller) error
        % estimate.
        errDecIdx = minErr > abs(Ctab{i-1, j-1});
        minErr(errDecIdx) = abs(Ctab{i-1, j-1}(errDecIdx));
        
        
        
        % If there is nothing to update we shall not continue the
        % iterations.
        if ~any(cGoodIdx(:))
            return;
        end
        
        % Update Jt estimate.
        updateIdx = cGoodIdx & errDecIdx;
        Jt(updateIdx) = Jtab{i,j}(updateIdx);
        
        % Update the power (p)
        p = p + 2;
    end
    
    % Update step length.
    h = h/2;
end
