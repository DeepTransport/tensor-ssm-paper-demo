function sams = priorsam(model, N)
sams = randn(model.d, N);
a = .4 + .6*normcdf(sams(1, :));
sams_x = randn(model.m, N);
sams = [sams; sams_x];
end