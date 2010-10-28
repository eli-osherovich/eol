function [x, allF, finalG, allGnorm, allX] = sesop_eo(x0, FuncAxStruct, FuncXStruct, options, params)


% maximal number of iterations
maxIter = options.maxSesopIter;

% maximal number of the Newton iterations
%maxNewtonIter = options.maxNewtonIter;

% number of previous gradients (including current one)
nPrevGrads = options.NSesopPrevGrads;

% number of previous steps
nPrevSteps = options.NSesopPrevSteps;

% misc directions (change it if you add code for misc. directions)
nMiscDirs = 0;

% problem dimensionality
N = numel(x0);


%% allocate space for previous gradients and directions along with those
% transformed by the lin. operators A_i

% previous gradients
PrevGrads = cell(1, nPrevGrads);
APrevGrads = cell(1, length(FuncAxStruct));

% previous steps
PrevSteps = cell(1, nPrevSteps);
APrevSteps = cell(1, length(FuncAxStruct));

% miscelaneous directions (see user code remark below)
MiscDirs = cell(1, nMiscDirs);
AMiscDirs = cell(1, length(FuncAxStruct));

% each linear operator has its own set of transformed directions
[APrevGrads{:}] = deal(cell(1, nPrevGrads));
[APrevSteps{:}] = deal(cell(1, nPrevSteps));
[AMiscDirs{:}] = deal(cell(1, nMiscDirs));


% union of the transformed directions
% (per linear operator)
ADirs = cell(1, length(FuncAxStruct));


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

x_newton = x0(:);
stepLength = 0;

%global NEWTGRAD;
%old_grad  = 0;

if ~isfield(options, 'display') || options.display
	fprintf('%10s %10s %15s %15s %15s\n','Iteration','FunEvals','Step Length','Function Val','Opt Cond');
end

% generate empty cell array of proper size (used by some functions)
empty = cell(size(FuncAxStruct));

for i = 0:maxIter
    %disp(i);
    
    % go to the new position
    x = x_newton;
    
    Ax = applyLinOpFwd(FuncAxStruct, x);
    
    % calculate function value and gradient
    [val, grad] = calc_EDx(x, Ax, FuncAxStruct, FuncXStruct, empty, []);
    gradNorm = norm(grad(:));
    
    %% use restart similar to the Polak–Ribiere method
    %     beta = grad(:)'*(grad(:)-old_grad(:));
    %     if beta <= eps
    %         restart('Negative beta');
    %     end
    %     old_grad = grad;
    
    % print statistics
	if ~isfield(options, 'display') || options.display
		fprintf('%10d %10d %15.5e %15.5e %15.5e\n',i,0,stepLength,val,gradNorm);
	end
    
    % record x, obj. function value, and gradient norm (if requested)
    if exist('allX', 'var'), allX(:,i+1) = x(:); end
    if isfield(options, 'dump') && options.dump
        save(fullfile(options.dumpdir, ['iter', num2str(i),'.mat']), 'x');
    end
    if exist('allF', 'var'), allF(i+1) = val; end
    if exist('allGnorm', 'var'), allGnorm(i+1) = norm(grad(:)); end
    
    % apply preconditioner to the gradient
    if   isfield(options, 'precond') && options.precond
        grad = options.mult_precond(grad,[], [], [], options.par);
    end

    
    % update previous gradients g, as well as A*g
    PrevGrads(1:end-1) = PrevGrads(2:end);
    PrevGrads{end} = -grad(:)/(gradNorm+eps);
    
    
    % calculate misc. dirs
    % place your code here
    %
    %
    % 
    % coarse grid direction
    if options.useCoarseGridDirection && i>0 && mod(i, options.useCoarseGridDirectionEvery) == 0
        cgd = GenerateCoarseGridDirection(x, grad, @RestrictNonCircularFullWeighting, @InterpolationNonCircularLinear, @CoarsenDeconvProblem, FuncAxStruct, FuncXStruct, options, params);
        MiscDirs = cgd(:)/(norm(cgd(:))+eps);
    else
        [AMiscDirs{:}] = deal(cell(1, nMiscDirs));
    end
    
    for k = 1:length(FuncAxStruct)
        APrevGrads{k}(1:end-1) = APrevGrads{k}(2:end);
        APrevGrads{k}{end} = applyLinOpFwd(FuncAxStruct, PrevGrads{end},k);
        AMiscDirs{k} = applyLinOpFwd(FuncAxStruct, MiscDirs,k);
    end
    
    % construct subspace
    Dirs = [PrevGrads, PrevSteps, MiscDirs];
    
    for k = 1:length(FuncAxStruct)
        ADirs{k} = [APrevGrads{k}, APrevSteps{k}, AMiscDirs{k}];
    end
    
    % calculate new position
    x_newton = newton_eo(x,Ax,FuncAxStruct,FuncXStruct, ADirs, Dirs, options);
    %[x_newton, allF_newton, finalG_newton, allGnorm_newton] = newtoneo(x,Ax,FuncAxStruct,FuncXStruct, ADirs, Dirs, options);
    %NEWTGRAD = horzcat(NEWTGRAD,allGnorm_newton);
    
    stepLength = norm(x_newton(:) - x(:));

    %     if stepLength/gradNorm < 0.01
    %            restart('short step');
    %     elseif nPrevSteps > 0
    if nPrevSteps > 0
        % update previous directions if we use them
        new_dir = x_newton - x;
        % normalize new direction
        new_dir = new_dir/(stepLength + eps);
        
        
        PrevSteps(1:end-1) = PrevSteps(2:end);
        PrevSteps{end} = new_dir;
        for k = 1:length(FuncAxStruct)
            APrevSteps{k}(1:end-1) = APrevSteps{k}(2:end);
            APrevSteps{k}{end} = applyLinOpFwd(FuncAxStruct, PrevSteps{end}, k);
        end
    end
end


if exist('finalG', 'var'), finalG=grad; end    
    
%     function [] = restart(msg)
%         disp(strcat(msg,':', 'resetting'))
%         PrevSteps(:) = 0;
%         PrevGrads(:) = 0;
%         for inner_idx = 1:length(FuncAxStruct)
%             APrevSteps{k}(:) = 0;
%             APrevGrads{k}(:) = 0;
%         end
%     end
end
