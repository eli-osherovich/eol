function t = step_eo(i, data, i_lo, i_hi, tMin, tMax)
% interpolation 


% Copyright 2010 Eli Osherovich.

% t0 is always t_lo
t0 = data(i_lo, 1);
f0 = data(i_lo, 2);
d0 = data(i_lo, 3);

t1 = data(i, 1);
f1 = data(i, 2);
d1 = data(i, 3);

% Pre-calculate cubic interpolation as it will be considered in most all
% cases.
tc = cubic_interpolate(t0, t1, f0, f1, d0, d1);

% Pre-calculate quadratic (secant) interpolation. It will be considered in
% certain cases only.
ts = t0 - d0*(t1 - t0)/(d1 - d0);

% First case. A higher function value. The minimum is bracketed. If the
% cubic step is closer to t_lo tan the quadratic step, the cubic step is
% taken, else, the average of the cubic and quadratic steps is taken.
% Quadratic step is well-defined in this case, i.e., the parabola is convex
% and has exactly one miminimum.
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

% Safe-guard the step.
% t = min(t, tMax);
% t = max(t, tMin);

% If bracket is found make sure it is not too close to its endpoints.
% if ~isempty(i_hi)
%     t2 = data(i_hi, 1);
%     if t2 > t0
%         t = min(t0 + 0.66*(t2-t0), t);
%     else
%         t = max(t0 + 0.66*(t2-t0), t);
%     end
% end
