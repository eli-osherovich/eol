function self = ctranspose(self)
    self.AdjointFlag = ~self.AdjointFlag;
    [self.RangeNumel, self.ImageNumel] = deal(self.ImageNumel, self.RangeNumel);
end
