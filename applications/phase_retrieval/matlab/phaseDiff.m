function pd = phaseDiff(phase1, phase2)
%PHASEDIFF - find phase difference

% the difference is normalized to the range [-pi, pi]
pd = phase1 - phase2;
pd = mod(pd + pi, 2*pi) - pi;
