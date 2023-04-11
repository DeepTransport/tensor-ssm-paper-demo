function y = tg_cdf(x, t)
%{
the truncated gaussian cdf for each element in x
[-t, t] to [0, 1]
-t gives value 0
t gives value 1
%}
if isscalar(t) && t>0
    C = 1/(1 - 2*normcdf(-t));
else
    error('tail should be positive scalar')
end
y = C * (normcdf(x) - normcdf(-t));
y(x==t) = 1;
y(x==-t) = 0;
y(x > t) = 1;
y(x < -t) = 0;
end