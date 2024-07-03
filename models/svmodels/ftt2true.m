function [gamma, tau, mu, phi, a, delta, nu1, nu2] = ftt2true(X, model)
%Transform the true parameters (potentially bounded) to ftt paras (unbounded)
% X is of size d * N
% outputs are 1*N vectors
N = size(X, 2);
gamma = normcdf(X(1,:));
tau = exp(X(2, :));

ind = 3;

% compute mu
if isnan(model.pre.para(3))
    mu = X(ind, :).*tau;
    ind = ind + 1;
else
    mu = repmat(model.pre.para(3), 1, N);
end

% compute phi
if isnan(model.pre.para(4))
    phi = X(ind, :).*tau;
    ind = ind + 1;
else
    phi = repmat(model.pre.para(4), 1, N);
end

% compute a
if isnan(model.pre.para(5))
    a = X(ind, :).*tau;
    ind = ind + 1;
else
    a = repmat(model.pre.para(5), 1, N);
end

% compute delta
if isnan(model.pre.para(6))
    delta = 4*(normcdf(X(ind, :)) - 0.5);
    ind = ind + 1;
else
    delta = repmat(model.pre.para(6), 1, N);
end

% compute nu1
if isnan(model.pre.para(7))
    nu1 = X(ind, :)*5 + 20;
    ind = ind + 1;
else
    nu1 = repmat(model.pre.para(7), 1, N);
end

% compute nu2
if isnan(model.pre.para(8))
    nu2 = X(ind, :)*5 + 20;
else
    nu2 = repmat(model.pre.para(8), 1, N);
end

end

