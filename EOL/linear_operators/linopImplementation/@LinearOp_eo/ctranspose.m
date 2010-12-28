function self = ctranspose(self)



% Copyright 2010 Eli Osherovich.



self.AdjointFlag = ~self.AdjointFlag;
[self.RangeNumelCurrent, self.ImageNumelCurrent] = ...
    deal(self.ImageNumelCurrent, self.RangeNumelCurrent);
end
