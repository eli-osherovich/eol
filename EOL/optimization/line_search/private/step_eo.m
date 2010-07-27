function t = step_eo(idx, data, idx_lo, idx_hi, tMin, tMax)
% Compute step-length. 



% Reference: 
%--------------------
% 1. Jorge J. Moré and David J. Thuente, “Line search algorithms
% with guaranteed sufficient decrease,” ACM Trans. Math. Softw. 20, no. 3
% (1994): 286-307.
%



% Copyright 2010 Eli Osherovich.

% t0 is always t_lo
t0 = data(idx_lo, 1);
f0 = data(idx_lo, 2);
d0 = data(idx_lo, 3);

t1 = data(idx, 1);
f1 = data(idx, 2);
d1 = data(idx, 3);

% Calculate cubic interpolation as it will be considered in most cases.
tc = cubic_interpolate(t0, t1, f0, f1, d0, d1);

% Calculate quadratic (secant) interpolation. It will be considered in
% certain cases only.
ts = t0 - d0*(t1 - t0)/(d1 - d0);

% First case. A higher function value. The minimum is bracketed. If the
% cubic step is closer to t_lo than the quadratic step, the cubic step is
% taken, else, the average of the cubic and quadratic steps is taken.
% Quadratic step is well-defined in this case, i.e., the parabola is convex
% and has exactly one minimum.
if f1 > f0
    tq = quadratic_interpolate(t0, t1, f0, f1, d0);
    
    if abs(tc - t0) < abs(tq - t0)
        t = tc;
    else
        t = tc + (tq - tc)/2;
    end
    
    % Second case. A lower function value and derivatives of opposite sign.
    % The minimum is bracketed. If the cubic step is closer to t_lo than
    % the quadratic (secant) step, the cubic step is taken, else the
    % quadratic step is taken.
elseif d1 * d0 < 0
    if abs(tc - t1) < abs(ts - t1)
        t = tc;
    else
        t = ts;
    end
    
    % Third case. A lower function value, derivatives of the same sign, and
    % the magnitude of the derivative decreases. The cubic step is only
    % used if the minimum of the cubic is beyond t1, i.e., the triplet
    % t_lo, t1, tc is ordered. The quadratic (secant)
    % step is also considered the step closest to t_lo is taken.
elseif abs(d1) < abs(d0)
    % unset tc if it does not extrapolate [t0, t1]
    if (tc - t0)*(tc - t1) < 0 
        tc = [];
    end
    if t1 > t0
        t = min([tc, ts, tMax]);
    else
        t = max([tc, ts, tMin]);
    end
else
    if t1 > t0
        t = tMax;
    else
        t = tMin;
    end
end


% If bracket is found make sure the step-length t is not too close to the
% endpoints of the bracket.
% The bracket should be given by [tMin, tMax] but we do not assume it and
% reconstruct it from the values t_lo and t_hi.
% if ~isempty(idx_hi)
%     t2 = data(idx_hi, 1);
%     tBrcktMin = min(t0, t2);
%     tBrcktMax = max(t0, t2);
%     % not too close from tBrcktMax
%     t = min(tBrcktMin + 0.9*(tBrcktMax - tBrcktMin), t);
%     % not too close from tBrcktMin
%     t = max(tBrcktMin + 0.1*(tBrcktMax - tBrcktMin), t);
% end
