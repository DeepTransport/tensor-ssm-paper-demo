function y = tg_inv(x, t)
%{
the truncated inverse gaussian cdf function for each element in x
[0, 1] to [-t, t]
0 gives value -t
1 gives value t
%}
if isscalar(t) && t>0
    C = 1/(1 - 2*normcdf(-t));
else
    error('tail should be positive scalar')
end
% temp = C\x + normcdf(-t);
y = norminv(C\x + normcdf(-t));
y(x<=0) = -t;
y(x>=1) = t;
end