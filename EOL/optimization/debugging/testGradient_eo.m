function [absErr, relErr] = testGradient_eo(x, FuncAxStruct, ...
    funcX, complexVarsFlag)
% TESTGRADIENT - test gradient (numerical vs. analytical).



% Copyright 2007-2010 Eli Osherovich.



% Calculate analytical gradient at the given point x.
Ax = applyMapping(FuncAxStruct, x);
[~, AnalGrad] = calcObjFunc(x, Ax, FuncAxStruct, funcX, complexVarsFlag);


% calculate numerical gradient 
NumGrad =  calcNumJacobian_eo(x, ...
    @(x) funcWrapper(x, FuncAxStruct, funcX, complexVarsFlag));


% Display results if no output requested (interactive mode).
if 0 == nargout
    figure;
    if complexVarsFlag
        subplot(411);plot(real([AnalGrad NumGrad])); title('Re(grad)')
        subplot(412);plot(imag([AnalGrad NumGrad])); title('Im(grad)');
        
        subplot(413);plot(real(AnalGrad - NumGrad)); title('Re(grad diff)');
        subplot(414);plot(imag(AnalGrad - NumGrad)); title('Im(grad diff)');
    else
        subplot(211);plot([AnalGrad NumGrad]); title('grad');
        subplot(212);plot(AnalGrad - NumGrad); title('grad diff');
    end
    fprintf('Max abs. error %e\n', max(abs(AnalGrad - NumGrad)));
    fprintf('Max rel. error %e\n', max(abs((AnalGrad -NumGrad)./AnalGrad)));
else
    absErr = max(abs(AnalGrad - NumGrad));
    relErr = max(abs((AnalGrad -NumGrad)./AnalGrad));
end

function val = funcWrapper(x, FuncAxStruct, funcX, complexVarsFlag)
Ax =  applyMapping(FuncAxStruct, x);
val = calcObjFunc(x, Ax, FuncAxStruct, funcX, complexVarsFlag);
