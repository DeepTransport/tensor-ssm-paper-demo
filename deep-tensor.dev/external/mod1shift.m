function Ps = mod1shift(P, delta)

delta = reshape(delta, size(P, 1), 1);
delta = repmat(delta, 1, size(P, 2));
Ps = mod(P+delta, 1);
