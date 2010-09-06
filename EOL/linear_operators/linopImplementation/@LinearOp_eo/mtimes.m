function Ax = mtimes(self, x)
    % If x is a linear opertor, the result is a linear operator
    % chain (which is, of course, a linear operator too).
    if isa(x, 'LinearOp_eo')
        Ax = LinearOpChain_eo({self, x});
        return;
    end
    
    % Check that the X's size is consistent with the rangeNumel of
    % the operator.
    if self.RangeNumel ~= numel(x)
        error('EOL:LinearOp:mtimes:wrongDimension', ...
            'Matrix dimensions must agree.');
    end
    
    % Apply operator in an appropriate mode, either forward or
    % adjoint.
    if self.AdjointFlag
        Ax = ApplyAdjoint(self, x);
    else
        Ax = ApplyForward(self, x);
    end
    
    if self.MinusFlag
        Ax = -Ax;
    end
end
