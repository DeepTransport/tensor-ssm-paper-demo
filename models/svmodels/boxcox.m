function [y, dydx] = boxcox(x, delta)
% compute the boxcox transformation of x given delta
% x and delta are 1*n vectors
% output y and dydx are also 1*n vectors
% dydx contains the derivative of x given delta
if length(x) ~= length(delta)
    error("X and delta should have same number of elements.")
end

ind1 = delta == 0;
ind2 = delta ~= 0;

y(ind1) = log(x(ind1));
y(ind2) = (x(ind2).^delta(ind2) - 1)./delta(ind2);

if nargout > 1
    dydx(ind1) = 1./x(ind1);
    dydx(ind2) = x(ind2).^(delta(ind2) -1);
end
end