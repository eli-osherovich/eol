function [phaseLo, phaseHi, phaseTrue] = genPhaseLimits_eo(x, phaseUncert,...
    complexVarsFlag)
%GENPHASELIMITS - generate boundaries of a phase uncertainty
% interval.
%
% [phaseLo, phaseHi] = genPhaseLimits_eo(x, phaseUncert) generate
% the lower (PHASELO) and upper (PHASEHI) bound of the phase
% uncertainty interval. The width (PHASEHI - PHASELO) of the
% interval equals to PHASE_UNCERT. In addition, it is assumed that x is a
% Fourier transform of a real signal, hence, it has certain symmetry
% (Hermitian) properties.
%
% [phaseLo, phaseHi] = genPhaseLimits_eo(x, phaseUncert, false) do not
% assume that the original signal is real.
%
% [phaseLo, phaseHi, phaseTrue] = genPhaseLimits_eo(...) returns
% the true phase as well.
%
% Implementation details:
%------------------------
% The interval bounds are not limited/normalized to any particular
% interval. It is guaranteed that PHASEHI-PHASELO=PHASE_UNCERT. The
% true phase PHASETRUE lies in the interval [-pi, pi] as returned
% by MATLAB's angle() function.
%
% The true phase is uniformly distributed in the interval 
% [PHASELO, PHASEHI].
%
% Shape of the output array is equal to that of X.
%
% PHASEUNCERT must be either a scalar (same uncertainty
% everywhere) or an array of shape equal to that of X.


% Copyright 2010 Eli Osherovich.
    
    
% Get true phase.
phaseTrue = angle(x);

% Generate lower bound.
phaseLo = phaseTrue - phaseUncert.*rand(size(phaseTrue));


% Generate upper bound (without the assumption that x is a Fourier
% transform of a real signal).
phaseHi = phaseLo + phaseUncert;

% Update the limits if the symmetry of a real origin is requested
if nargin < 3 || ~complexVarsFlag
    % in case of a real input, the Fourier transform is determined by
    % (approximately) half of the output values.
    % Note that the shape of X must be correct in this case.
    xSize = size(x);
    idx1 = 1:floor(numel(x)/2) + 1;
    subsCell = cell(1, ndims(x));
    [subsCell{:}] = ind2sub(xSize, idx1);
    for d = 1:ndims(x)
        subsCell{d} = xSize(d) - subsCell{d} + 2;
        subsCell{d}(subsCell{d} > xSize(d)) = 1;
    end
    idx2 = sub2ind(xSize, subsCell{:});
    
    % Some values are ought to be real. We shall not change the boundaries
    % for these.
    nonrealIdx = idx1 ~= idx2;
    
    phaseLoTmp = phaseLo;
    phaseHiTmp = phaseHi;
    phaseLo(idx2(nonrealIdx)) = -phaseHiTmp(idx1(nonrealIdx));
    phaseHi(idx2(nonrealIdx)) = -phaseLoTmp(idx1(nonrealIdx));
end
