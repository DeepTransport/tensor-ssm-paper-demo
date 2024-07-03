function [gamma, tau, mu, phi, a, delta, nu1, nu2] = true2ftt(X, model)
%Transform the true parameters (potentially bounded) to ftt paras (unbounded)
% X is of size d * N
% outputs are 1*N vectors
N = size(X, 2);
gamma = norminv(X(1,:));
tau = log(X(2, :));
ind = 3;

% compute mu
if isnan(model.pre.para(3))
    mu = X(ind, :)./X(2, :);
    ind = ind + 1;
else
    mu = model.pre.para(3)./X(2, :);
end

% compute phi
if isnan(model.pre.para(4))
    phi = X(ind, :)./X(2, :);
    ind = ind + 1;
else
    phi = model.pre.para(4)./X(2, :);
end

% compute a
if isnan(model.pre.para(5))
    a = X(ind, :)./X(2, :);
    ind = ind + 1;
else
    a = model.pre.para(5)./X(2, :);
end

% compute delta
if isnan(model.pre.para(6))
    delta = norminv(X(ind, :)/4 + 0.5);
    ind = ind + 1;
else
    delta = repmat(norminv(model.pre.para(6)/4 + 0.5), 1, N);
end

% compute nu1
if isnan(model.pre.para(7))
    nu1 = X(ind, :)/5 - 4;
    ind = ind + 1;
else
    nu1 = repmat(model.pre.para(7)/5 - 4, 1, N);
end

% compute nu2
if isnan(model.pre.para(8))
    nu2 = X(ind, :)/5 - 4;
else
    nu2 = repmat(model.pre.para(8)/5 - 4, 1, N);
end

end

