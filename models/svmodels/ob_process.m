function Y = ob_process(model, thetax, ~)
% this is only used to generate observations so use true theta
% X only contains states
% gamma = 2 * normcdf(thetax(1, :)) - 1;
% mu = thetax(2, :);
% tau = exp(thetax(3, :));
% phi = thetax(4, :);
[gamma, tau, mu, phi, a, delta, nu1, nu2] = ftt2true(thetax(1:model.d, :), model);

normstd = sqrt(boxcoxinv(thetax(model.d+1,:).*tau, delta));
Y = normstd.* trnd(nu2, 1, size(thetax,2));

end