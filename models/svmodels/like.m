function pdf = like(model, thetax, t)
[gamma, tau, mu, phi, a, delta, nu1, nu2] = ftt2true(thetax(1:model.d, :), model);

y = model.Y(:, t);
normstd = sqrt(boxcoxinv(thetax(model.d+1,:).*tau, delta));

pdf = tpdf(y./normstd, nu2)./normstd;
pdf(isnan(pdf)) = 0;
pdf(nu2<4) = 0;


end