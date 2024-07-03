function sams = priorsam(model, N)
sams12 = zeros(2, N);
gamma = 2 * betarnd(20, 1.5, 1, N) - 1;
sams12(1,:) = norminv(max(min(gamma, 1-1e-10), 1e-10));
% tau2 = random('InverseGaussian', 1, 0.005, [1, N]);
tau2 = 1./ gamrnd(1, 1/0.005, 1, N);
sams12(2,:) = log(sqrt(tau2));


if isnan(model.pre.para(3))
    mu = normrnd(0, sqrt(5), [1, N]);
    sams3 = mu;
    ind = 1;
else
    mu = model.pre.true(3)./ sqrt(tau2);
    sams3 = [];
    ind = 0;
end

if isnan(model.pre.para(4)) && isnan(model.pre.para(5))
    sams45 = normrnd(0, sqrt(2), [2, N]);
    ind = ind + 2;
elseif isnan(model.pre.para(4)) || isnan(model.pre.para(5))
    sams45 = normrnd(0, sqrt(2), [1, N]);
    ind = ind + 1;
else
    sams45 = [];
    ind = ind;
end

ind678 = model.d - 2 - ind;
if ind678 > 0
    sams678 = normrnd(0, 1, [ind678, N]);
else
    sams678 = [];
end


samsx = normrnd(mu, sqrt(1./(1-gamma.^2)));
sams = [sams12;sams3;sams45;sams678;samsx];

% sams = randn(model.d+model.m, N);
% sams(end,:) = sams(end, :) + model.theta(2);

end