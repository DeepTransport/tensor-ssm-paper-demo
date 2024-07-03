% test ftt2true and true2 ftt

model.pre.para = [NaN, NaN, NaN, -2, NaN, -5, NaN, 7];
X = randn(5, 10);
[gamma, mu, tau, phi, a, delta, nu1, nu2] = ftt2true(X, model);
[gammax, mux, taux, phix, ax, deltax, nu1x, nu2x] = true2ftt([gamma; mu; tau; a; nu1], model);







