function d = lbfgsDir_eo(validIdx, wrapAround, grad, PrevSteps, PrevGrads, H0)

nPrev = length(PrevSteps);



% number of valid steps in the memory
if wrapAround
    nValid = nPrev;
    bwdIdx = [validIdx:-1:1, nPrev:-1:validIdx+1];
    fwdIdx = [validIdx+1:nPrev, 1:validIdx];
else
    nValid = validIdx;
    bwdIdx = validIdx:-1:1;
    fwdIdx = 1:validIdx;
end

% preallocate space
rho = zeros(nValid, 1);
alpha = zeros(nValid, 1);
beta = zeros(nValid, 1);


d = grad;


for i = bwdIdx
    rho(i) = 1/(PrevSteps{i}'*PrevGrads{i});
    alpha(i) = (PrevSteps{i}'*d) * rho(i);
    d = d - alpha(i) * PrevGrads{i};
end

% mutiply by the initial Hessian approximation (scalar matrix in our case).
d = H0 * d;

for i = fwdIdx
    beta(i) = (PrevGrads{i}' * d) * rho(i);
    d = d + PrevSteps{i} * (alpha(i)-beta(i));
end

