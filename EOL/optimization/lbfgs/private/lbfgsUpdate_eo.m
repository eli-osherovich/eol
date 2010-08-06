function [PrevSteps, PrevGrads,H0, validIdx, wrapAround] = lbfgsUpdate_eo(validIdx, wrapAround,...
    x_diff, grad_diff, PrevSteps, PrevGrads,H0)
% LBFGSUPDATE - update the LBFGS memory.
% 
% The function updates the LBFGS memory, i.e., stores the most recent step
% and gradient difference. Additionally it updates the scalar (matrix) H0
% which is used as an initial approximation to the inverse Hessian matrix.
%
%
% Reference:
% ----------
% Numerical Optimization by J. Nocedal and S.J. Wright.



% Copyright 2010 Eli Osherovich.


gradstep = real(grad_diff'*x_diff);

if gradstep > 1e-10
    % Update is done only the above product is bounded from zero.
    % Otherwise the inverse Hessian may not be positive definite.
    nPrev = length(PrevSteps);

    if validIdx >= nPrev % it should never be greater than...
        wrapAround = 1;
        validIdx = 1;
    else
        validIdx = validIdx + 1;
    end
    PrevSteps{validIdx} = x_diff;
    PrevGrads{validIdx} = grad_diff;

    % Update the inverse Hessian approximation (a scalar).
    H0 = gradstep/real(grad_diff'*grad_diff);
end
