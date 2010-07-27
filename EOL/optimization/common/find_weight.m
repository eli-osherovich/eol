function w = find_weight(F)
if ~isfield(F, 'weight') || isempty(F.weight)
    w = 1;
else
    w = F.weight;
end
