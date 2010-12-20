function Ax = mtimes(self, x)
    % If x is a linear opertor, the result is a linear operator
    % chain (which is, of course, a linear operator too).
    if isa(x, 'LinearOp_eo')
        Ax = LinearOpChain_eo({self, x});
    elseif isnumeric(x)
        Ax = applyOp(self, x);
    elseif iscell(x)
        % pre-allcote space
        Ax = cell(size(x));
    
        for i = 1:numel(x)
            % Skip if empty
            if isempty(x{i})
                continue;
            end
            Ax{i} = applyOp(self, x);
        end
    else
        error('EOL:PenaltyFunc:WrongArgType', ...
            ['x must be a linear operator, numeric vector, or cell array.',...
            '\nInstead got %s.'], class(x));
    end

    function Ax = applyOp(op, x)
        % Check that the X's size is consistent with the rangeNumel of
        % the operator.
        if op.RangeNumel ~= numel(x)
            error('EOL:LinearOp:mtimes:wrongDimension', ...
                'Matrix dimensions must agree.');
        end
        
        % Apply operator in an appropriate mode, either forward or
        % adjoint.
        if op.AdjointFlag
            Ax = ApplyAdjoint(op, x);
        else
            Ax = ApplyForward(op, x);
        end
        
        if op.MinusFlag
            Ax = -Ax;
        end
    end
end
