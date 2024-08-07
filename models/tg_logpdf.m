function y = tg_logpdf(x, t)
%{
the multivariate truncated gaussian pdf for each column in x
[-t, t] to [0,some value]
-t t should give values close to 0
%}
if isscalar(t) && t>0
    C = 1/(1 - 2*normcdf(-t));
else
    error('tail should be positive scalar')
end
y = log(C) + (-0.5 * x.^2) - log(sqrt(2*pi));
% y(x > t) = 0;
% y(x < -t) = 0;
y = sum(y, 1);

end