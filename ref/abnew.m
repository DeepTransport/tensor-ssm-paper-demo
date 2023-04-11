function [L, c] = abnew(sams, sdnum)
% flag chooses different reference dis
% quant is the quantile, sams are samples of size d*N

n = size(sams, 2);
d = size(sams, 1);
mu = zeros(d, 1);
L = zeros(d);
for k = 1 : d
    data = sams(k, :)';
    sd = sqrt(var(data));
    mu = mean(data);
    x1 = mu - sdnum * sd;
    x2 = mu + sdnum * sd;
    L(k, k) = x2-x1;
    c(k) = mu;
end
c=c';
end