function d = lbfgsDirPersist_eo(validIdx, wrapAround, grad, PrevSteps, PrevGrads, H0)

% NOTE!
% due to pesistent memory, only one instance of LBFGS may run at
% any given moment.


persistent rho nPrev_old



nPrev = length(PrevSteps);


% size change or first time use
if isempty(nPrev_old) || nPrev_old ~= nPrev
    rho = zeros(1, nPrev);
    nPrev_old = nPrev;
end

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
alpha = zeros(nValid, 1);
beta = zeros(nValid, 1);


d = grad;

if 0 ~= validIdx
    rho(validIdx) = 1/real(PrevSteps{validIdx}'*PrevGrads{validIdx});
end

for i = bwdIdx
    alpha(i) = real(PrevSteps{i}'*d) * rho(i);
    d = d - alpha(i) * PrevGrads{i};
end

% mutiply by the initial Hessian approximation (scalar matrix in our case).
d = H0 * d;

for i = fwdIdx
    beta(i) = real(PrevGrads{i}'*d) * rho(i);
    d = d + PrevSteps{i} * (alpha(i)-beta(i));
end
