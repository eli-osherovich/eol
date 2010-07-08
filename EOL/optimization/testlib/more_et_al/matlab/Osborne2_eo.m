function [fvec, Jac, hessV, hessVfunc, x0_std] = Osborne2_eo(x, ~)
%OSBORNE2 - Osborne 2 function.

% Number of variables: N = 11
% Number of equations: M = 65
% Standard starting point x0_std = [1.3; 0.65; 0.65; 0.7; 0.6; 3;
% 5; 7; 2; 4.5; 5.5]
% Min function value = 4.01377E-2
% Argmin x =


% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [1.3; 0.65; 0.65; 0.7; 0.6; 3; 5; 7; 2; 4.5; 5.5];
M = 65;
N = 11;

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end

% initialize some constants
t = ((1:M)' - 1)/10;
y = [1.366; 1.191; 1.112; 1.013; 0.991; 0.885; 0.831; 0.847; 0.786; ...
     0.725; 0.746; 0.679; 0.608; 0.655; 0.616; 0.606; 0.602; 0.626; ...
     0.651; 0.724; 0.649; 0.649; 0.694; 0.644; 0.624; 0.661; 0.612; ...
     0.558; 0.533; 0.495; 0.500; 0.423; 0.395; 0.375; 0.372; 0.391; ...
     0.396; 0.405; 0.428; 0.429; 0.523; 0.562; 0.607; 0.653; 0.672; ...
     0.708; 0.633; 0.668; 0.645; 0.632; 0.591; 0.559; 0.597; 0.625; ...
     0.739; 0.710; 0.729; 0.720; 0.636; 0.581; 0.428; 0.292; 0.162; ...
     0.098; 0.054];


fvec = y - (x(1)*exp(-t*x(5)) + x(2)*exp(-(t-x(9)).^2*x(6)) ...
            + x(3)*exp(-(t-x(10)).^2*x(7)) + x(4)*exp(-(t-x(11)).^2*x(8)));



if nargout > 1 % Jacobian requested
    Jac = zeros(M, N);
    Jac(:, 1) = -exp(-t*x(5));
    Jac(:, 2) = -exp(-(t-x(9)).^2*x(6));
    Jac(:, 3) = -exp(-(t-x(10)).^2*x(7));
    Jac(:, 4) = -exp(-(t-x(11)).^2*x(8));
    Jac(:, 5) = x(1)*t.*exp(-t*x(5));
    Jac(:, 6) = x(2)*(t-x(9)).^2.*exp(-(t-x(9)).^2*x(6));
    Jac(:, 7) = x(3)*(t-x(10)).^2.*exp(-(t-x(10)).^2*x(7));
    Jac(:, 8) = x(4)*(t-x(11)).^2.*exp(-(t-x(11)).^2*x(8));
    Jac(:, 9) = -2*x(2)*x(6)*(t-x(9)).*exp(-(t-x(9)).^2*x(6));
    Jac(:, 10) = -2*x(3)*x(7)*(t-x(10)).*exp(-(t-x(10)).^2*x(7));
    Jac(:, 11) = -2*x(4)*x(8)*(t-x(11)).*exp(-(t-x(11)).^2*x(8));
    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:Osborne2:HessianNotSupported', ['Hessian or ' ...
                            'Hessian-vector product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
