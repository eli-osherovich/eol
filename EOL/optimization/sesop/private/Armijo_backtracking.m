function t = Armijo_backtracking(x0, Ax0, val0, grad0, func_Ax_Struct, func_x_Struct, Adir, dir, options)

% Reference: Convex Optimization by Stephen Boyd and Lieven Vandenberghe.
% Our notion follows that of the book.

% alpha - measure of "function decrease sufficiency".
% it usually lies in the interval (0, 0.5)
alpha = options.LineSearchDecreaseFactor;


% initial t value is one - to fit Newton method best
t = 1;

% maximimal iterations
maxIter = 25;


% inner product between the gradient at x0 and direction 'dir'
%ip = grad0(:)'*dir(:);

Ax = Ax0; % equivalent to t = 1

for i = 1:maxIter
    
    % next step trial
    x = x0 + t*dir;
    for k = 1:length(func_Ax_Struct)
        Ax{k} = Ax0{k} + t * Adir{k};
    end
    
    % uncomment this is gradien calucations are cheap
    %[val, grad] = calc_EDx(x, Ax, func_Ax_Struct, func_x_Struct, Adir, dir);
    val = calc_EDx(x, Ax, func_Ax_Struct, func_x_Struct, [], []);
    
    val0eps = eps*val0;
    
    if val > val0 + alpha*t*grad0 - val0eps - eps*val
        [val, grad] = calc_EDx(x, Ax, func_Ax_Struct, func_x_Struct, Adir, dir, true);
        t = cubic_interpolate(0, t, val0, val, grad0, grad);
        %t = quadratic_interpolate(0, t, val0, val, gradDotDir);
        %t = t/2;
    else
        return
    end
end
