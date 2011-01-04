function self = conj(self)


% Copyright 2010 Eli Osherovich.



self.AdjointFlag = ~self.AdjointFlag;
[self.RangeNumel, self.ImageNumel] = ...
    deal(self.ImageNumel, self.RangeNumel);
end
