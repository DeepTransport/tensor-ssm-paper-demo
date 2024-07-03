function [y, dydx] = boxcoxinv(x, delta)
% compute the boxcox transformation of x given delta
% x and delta are 1*n vectors
% output y and dydx are also 1*n vectors
% dydx contains the derivative of x given delta
if length(x) ~= length(delta)
    error("X and delta should have same number of elements.")
end

ind1 = delta == 0;
ind2 = delta ~= 0;

y(ind1) = exp(x(ind1));
y(ind2) = (1+delta(ind2).*x(ind2)).^(1./delta(ind2));

ind3 = (delta ~= 1) & (delta.*x < -1);
y(ind3) = NaN;

if nargout > 1
    dydx(ind1) = exp(x(ind1));
    dydx(ind2) = (1+delta(ind2).*x(ind2)).^(1./delta(ind2)-1);
    dydx(ind3) = NaN;
end
end