function [PrevSteps, PrevGrads,H0, validIdx, wrapAround] = lbfgsUpdate_eo(validIdx, wrapAround,...
    x_diff, grad_diff, PrevSteps, PrevGrads,H0)


gradstep = real(grad_diff'*x_diff);

if gradstep > 1e-10
    nPrev = length(PrevSteps);

    if validIdx >= nPrev % it should never be greater then...
        wrapAround = 1;
        validIdx = 1;
    else
        validIdx = validIdx + 1;
    end
    PrevSteps{validIdx} = x_diff;
    PrevGrads{validIdx} = grad_diff;

    % Update Hessian (scalar matrix) approximation.
    H0 = gradstep/real(grad_diff'*grad_diff);
end
