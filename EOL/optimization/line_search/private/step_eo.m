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
% efficiently make it infeasable at later stages.
if ~isreal(tc)
    tc = [];
end

if ~isempty(idx_hi)
    % Bracket is known. Use interpolation.
    
    % Unset tc it lies outside the bracket [t0, t1]
    if (tc - t0)*(tc - t1) > 0
        tc = [];
    end
    
    tq = quadratic_interpolate(t0, t1, f0, f1, d0);
    if (tq - t0)*(tq - t1) > 0
        tq = [];
    end
    
    % Some interior points.
    tThird = t0 + (t1 - t0)/3;
    tMid = t0 + (t1 - t0)/2;
    
    % First case: we have a "promising" bracket with derivatives of
    % different sign at its endpoints. Cubic interpolation will always give
    % a meaningful result in this case. Hence we use it.
    if d1 * d0 < 0;
        t = tc;
        
    elseif f1 > f0
        % Second case: the derivatives have the same sign and F1 > F0.
        % Quadratic interpolation cannot be correct here but we still
        % consider it.
        
        
        % Use the point which is farthers from T0.
        if abs(tc - t0) < abs(tq - t0)
            t = tq;
        else
            t = tc;
        end
        % But never go too far from T0.
        t = min(t, tMid);
        
    else
        % Use the point which is closest to T0.
        if abs(tc - t0) < abs(tq - t0)
            t = tc;
        else
            t = tq;
        end
        % But never go too far from T0.
        t = min(t, tThird);
    end
    
    % Safeguarding: do not let t be too close to the bracket endpoints.
    t = min(t, tMax - 0.1 * (tMax - tMin));
    t = max(t, tMin + 0.1 * (tMax - tMin));
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
        % Quadratic (secant) approximation.
        ts = t0 - d0*(t1 - t0)/(d1 - d0);
        % Bracket is not known T0 > T1, F0 < F1, and ABS(D0) <
        % ABS(D1).TC is considered only if it lies beyond T0. TS is also
        % considered since always lies on the correct side.
        t = min([ts, tc, tMax]);
    else
        t = min([tc, tMax]);
    end
end

