function pdfval = transition(model, thetax, t)
% nx = size(thetax,2);
if t == 1
    y = 0;
else
    y = model.Y(t-1);
end

[gamma, tau, mu, phi, a, delta, nu1, nu2] = ftt2true(thetax(1:model.d, :), model);

normmu =  mu + gamma.*(thetax(model.d+model.m+1, :)-mu) + ...
    (phi ./ sqrt(boxcoxinv(thetax(model.d+model.m+1,:).*tau, delta))) * y +...
     a.* y^2;

pdfval = tpdf(thetax(model.d+1,:) - normmu, nu1);
pdfval(isnan(pdfval)) = 0;
pdfval(nu1<4) = 0;
end