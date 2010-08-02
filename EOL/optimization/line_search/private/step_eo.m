function t = step_eo(data, idx, idx_lo, idx_hi, tMin, tMax)
% Compute next step-length. 
% 
% Inputs:
%------------
% DATA - an Nx3 array which contains all information about the
% one-dimensional function f(x0 + t*dir) obtained so far. Each row contains
% a triplet [t, f, d] where t is the step length, f - function value and d
% - directional derivative (grad(f)'*dir) value.
%
% IDX - index of the current (latest) attempted step-length.
%
% IDX_LO - index of the step-length that yielded the lowest function value
% that satisfies the sufficient decrease condition. There are may be other
% step-length with even lower function values that, however, do no satisfy
% the sufficient decrease condition.
%
% IDX_HI - in case the minimum of f(x0 + t*dir) was bracketed this variable
% contains the index of the second end-point (the first end-point is stored
% in IDX_LO).
%
% TMIN - minimum allowed step-length.
%
% TMAX - maximum allowed step-length.



% Reference: 
%--------------------
% 1. Jorge J. Moré and David J. Thuente, “Line search algorithms
% with guaranteed sufficient decrease,” ACM Trans. Math. Softw. 20, no. 3
% (1994): 286-307.
%
% With significant changes.



% Copyright 2010 Eli Osherovich.

% t0 is always t_lo
% t1 is chosen differently depending on whether the bracket is found or
% not.
t0 = data(idx_lo, 1);
f0 = data(idx_lo, 2);
d0 = data(idx_lo, 3);

if ~isempty(idx_hi)
    % Bracket is known: use its endpoints.
    t1 = data(idx_hi, 1);
    f1 = data(idx_hi, 2);
    d1 = data(idx_hi, 3);
else
    % Bracket is not known: use another point which is closest to t_lo.
    [~, idx_tmp] =min(abs(data([1:idx_lo-1, idx_lo+1:idx], 1) - ...
                          data(idx_lo,1)));
    t1 = data(idx_tmp, 1);
    f1 = data(idx_tmp, 2);
    d1 = data(idx_tmp, 3);
end


% Calculate cubic interpolation as it will be considered in most cases.
tc = cubic_interpolate(t0, t1, f0, f1, d0, d1);
% In certain cases the cubic polynomial may have no minimum. Such cases are
% idicated by tc being complex. Setting tc to an empty array will
% efficiently eliminate it if later stages.
if ~isreal(tc)
    tc = [];
end

% Calculate quadratic (secant) interpolation. It will be considered in
% certain cases only.
ts = t0 - d0*(t1 - t0)/(d1 - d0);

if ~isempty(idx_hi)
    % Bracket is known. Use interpolation.
    tq = quadratic_interpolate(t0, t1, f0, f1, d0);
    % Unset tc  and tq if they lie outside the bracket [t0, t1]
    if (tc - t0)*(tc - t1) > 0
        tc = [];
    end
    if (tq - t0)*(tq - t1) > 0
        tq = [];
    end
    
    % Calculate mid point as it will be used in some cases.
    tMid = t0 + (t1 - t0)/2;
    
    if f1 > f0
        % First case: F1 > F0. Use TC if it is closer to T0 than TQ.
        % Otherwise use the average of TC and TQ.
        if abs(tc - t0) < abs(tq - t0)
            t = tc;
        else
            t = tc + (tq - tc)/2;
        end
    elseif d1 * d0 < 0
        % Second case: F1 <= F0 and the derivatives D0 and D1 have opposite
        % sign. Use TC if it is closer to T1. Otherwise use TS which is
        % guaranteed to lie within the bracket [T0, T1].
        if abs(tc - t1) < abs(ts - t1)
            t = tc;
        else
            t = ts;
        end
    else
        % Third case: F1 <= F0, the derivatives D0 and D1 have the same
        % sign. Use TC or the midpoint TMID (whichever lies closer to t0).
        % Proximity to T0 is preffered because T1 cannot satisfy the
        % sufficient decrease condition.
        if t1 > t0
            t = max([tc, tMid]);
        else
            t = min([tc, tMid]);
        end
    end
    t = min(t, tMax);
    t = max(t, tMin);
else
    % Bracket is not known. Use extrapolation.
    % Due to how WOLFELS_EO works this situation is only possible when each
    % tried T up to now resulted in a lower F value all of which satisfy
    % the sufficient decrease condition. Hence, we assume that all T's in
    % the DATA array constitute a monotonically increasing sequence, while
    % F's are monotonically decreasing.
    
    % Unset tc if it lies inside [t0, t1]
    if tc < t0
        tc = [];
    end
    if abs(d0) < abs(d1)
        % Fourth case: bracket is not known T0 > T1, F0 < F1, and ABS(D0) <
        % ABS(D1).TC is considered only if it lies beyond T0. TS is also
        % considered since always lies on the correct side.
        t = min([ts, tc, tMax]);
    else
        t = min([tc, tMax]);
    end
    % Do not let T be too close to T0.
    t = max(tMin, t);
end

