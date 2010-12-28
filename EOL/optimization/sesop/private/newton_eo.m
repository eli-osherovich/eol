function [x, allF, finalG, allGnorm, allX] = newton_eo(x0, Ax0, ...
        FuncAxStruct, funcX, AD, D, options)

% maximal number of iterations
maxIter = options.maxNewtonIter;

% scale factor of small Hessian eigenvalues
hessCond = options.HessianConditionNumber;

% gradient norm termination criteria
gradNormTerm = options.NewtonGradNormTermination;

% step size termination
stepSizeTerm = options.NewtonStepSizeTermination;

% directional derivative termination
directionalDerivTerm = options.directionalDerivTerm;

% gradient decrease termination
gradientDecreaseRatio = options.gradientDecreaseRatio;

% recalculate true Hessian every
recalcHess = options.recalcHessian;


% debug flag
debug = options.debug;

% problem dimensionality
N = numel(x0);


Adir = cell(1, length(FuncAxStruct));

% allocate space for output requested
if nargout > 1,
    % f(x_i) values. i-th column contains f_i
    allF = NaN(1, maxIter+1);
    if nargout > 2
        % return final gradient
        finalG = NaN(size(x0));
        if nargout > 3
            % norm(grad(x_i)) values. i-th column contains ||g_i||
            allGnorm = NaN(1, maxIter+1);
            if nargout > 4
                % collection of all x values. i-column row contains x_i
                allX = NaN(N, maxIter+1);
            end
        end
    end
end


x_new = x0;
Ax = Ax0;



for i = 1:maxIter
    
    earlyExit = false;
    
    if i > 1
        s = t*dir_alpha(:);
    end
    x = x_new;
    
        
    % compute value, gradient, and Hessian
    if (i == 1) || (mod(i,recalcHess) == 0)
        if (debug)
            fprintf('Re-calculating true Hessian\n');
        end
        [val, grad_alpha, hess_alpha] = calc_EDx(x, Ax, FuncAxStruct, funcX, AD, D, true);
        
        if (i == 1), % remember intial gradient
            grad_alpha0norm = norm(grad_alpha(:));
        end
        exactHessianFlag  = true;
    else
        exactHessianFlag = false;
        old_grad_alpha = grad_alpha;
        [val, grad_alpha] = calc_EDx(x, Ax, FuncAxStruct, funcX, AD, D, true);
        y = grad_alpha(:) - old_grad_alpha(:);
        
        sy = s'*y;
        if sy > 1e-10, % this test is not necessary once we have Wolfe linesearch
            %rho = 1/(s'*y);
            %tmp = eye(numel(grad_alpha)) - rho*y*s';
            %hess_alpha_inv = tmp'*hess_alpha_inv*tmp + rho*s*s';
            Hy = hess_alpha_inv*y;
            hess_alpha_inv = hess_alpha_inv + (1 + y'*Hy/sy)*(s*s')/sy-...
                (s*Hy' + Hy*s')/sy;
            % DFT update
            %tmp = hess_alpha_inv*s;
            %hess_alpha_inv = hess_alpha_inv + y*y'*rho - tmp*tmp'/(s'*tmp);
        else
            % do nothing, use previous hessian
            if (debug)
                fprintf('Skipping Hessian update\n');
            end
        end
    end

    
    % record x, obj. function value, and gradient norm (if requested)
    if exist('allX', 'var'), allX(:,i) = x(:); end
    if exist('allF', 'var'), allF(i) = val; end
    if exist('allGnorm', 'var'), allGnorm(i) = norm(grad_alpha(:)); end
    
    % check if gradient is small enough to terminate
    if norm(grad_alpha(:))/grad_alpha0norm < gradientDecreaseRatio
        if (debug)
            fprintf('Exiting on small gradient ratio\n');
        end
        earlyExit = true;
        break;
    end
    if norm(grad_alpha(:)) < gradNormTerm
        if (debug)
            fprintf('Exiting on small gradient norm\n');
        end
        earlyExit = true;
        break;
    end
    
    if (exactHessianFlag)
        % hessian may become non-symmetric due to numerical errors
        % hess_alpha = real((hess_alpha + hess_alpha')/2);
    
        % eigen decomposition of the Hessian
        [evect,eval] = eig(hess_alpha);
        d = diag(eval);
        
        assert (isreal(d) && isreal(evect));
    
        
        % replace negative values with positive ones of the same magnitude
        d = abs(d);
        
        % force the Hessian to be *sufficiently* positive definite
        %d=max(d,hessCond*max(d));
        valid_idx = d > numel(grad_alpha(:))*eps(max(d));

        d_inv = zeros(size(d));
        d_inv(valid_idx) = 1./d(valid_idx);
        % calculate Newton direction
        hess_alpha_inv = evect*diag(d_inv)*evect';
        %hess_alpha_inv = pinv(hess_alpha);
        %dir_alpha =  -evect*((evect'*grad_alpha(:))./d);
        dir_alpha = -hess_alpha_inv*grad_alpha(:);
    else
        dir_alpha = -hess_alpha_inv*grad_alpha(:);
    end
    
    % translate direction from the subspace into the space of x and Ax
    dir = 0;
    non_empty_idx = find(~cellfun('isempty',D));
    for k = 1:numel(dir_alpha)
        dir = dir + D{non_empty_idx(k)}*dir_alpha(k);
    end
    
    % test directional derivative
    if abs(grad_alpha'*dir_alpha) < directionalDerivTerm,
        if (debug)
            fprintf('Exitting on small directional derivative\n');
        end
        earlyExit = true;
        break;
    end
    
    
    
    %dir = reshape(dir, size(x0));
    for n = 1:length(FuncAxStruct)
        non_empty_idx = find(~cellfun('isempty',AD{n}));
        adir_tmp = 0;
        for k = 1:numel(dir_alpha)
            adir_tmp = adir_tmp + AD{n}{non_empty_idx(k)}*dir_alpha(k);
        end
        Adir{n} = adir_tmp;
    end
    
    % calulate step length
    t = Armijo_backtracking(x, Ax, val, grad_alpha(:)'*dir_alpha, FuncAxStruct, funcX, Adir, dir, options);
    %t = Wolfe_cubic(x, Ax, val, grad_alpha(:)'*dir_alpha, FuncAxStruct, funcX, Adir, dir, options);

    % terminate if step size is too small 
    if t*norm(dir(:)) < stepSizeTerm;
        if (debug)
            fprintf('Exitting on small step length\n');
        end
        earlyExit = true;
        break;
    end
        
    % perform step
    x_new = x + t*dir;
    for k = 1:length(FuncAxStruct)
        Ax{k} = Ax{k} + t*Adir{k};
    end
end

if ~earlyExit,    %update last step

    if debug
        fprintf('Max number of iterations reached\n');
    end

    switch nargout
        case 1
            x = x_new;
        case 2
            x = x_new;
            allF(i+1) = calc_EDx(x, Ax, FuncAxStruct, funcX, [], []);
        case 3
            x = x_new;
            [val, grad_alpha] = calc_EDx(x, Ax, FuncAxStruct, funcX, AD, D, true);
            allF(i+1) = val;
            finalG = grad_alpha;
        case 4
            x = x_new;
            [val, grad_alpha] = calc_EDx(x, Ax, FuncAxStruct, funcX, AD, D, true);
            allF(i+1) = val;
            finalG = grad_alpha;
            allGnorm(i) = norm(grad_alpha(:));
        case 5
            x = x_new;
            [val, grad_alpha] = calc_EDx(x, Ax, FuncAxStruct, funcX, AD, D, true);
            allF(i+1) = val;
            finalG = grad_alpha;
            allGnorm(i+1) = norm(grad_alpha(:));
            allX(i+1) = x(:);
        otherwise
            error('Oooops');
    end
else % exit by break
    if nargout > 2
        finalG = grad_alpha;
    end
end


