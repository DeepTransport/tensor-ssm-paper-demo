function X_new = st_process(model, thetax, t)
nx = size(thetax, 2);
X_new = zeros(1, nx);
[gamma, tau, mu, phi, a, delta, nu1, nu2] = ftt2true(thetax(1:model.d, :), model);

if t == 1
    y = 0;
else
    y = model.Y(t-1);
end


X_new(1, :) =  mu + gamma.*(thetax(model.d+1, :)-mu) + ...
    (phi ./ sqrt(boxcoxinv(thetax(model.d+1,:).*tau, delta))) * y +...
     a.* y^2 + trnd(nu1, 1, nx);
X_new(isnan(X_new)) = inf;
end