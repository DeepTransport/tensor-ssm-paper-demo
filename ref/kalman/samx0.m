function sams = samx0(thetas, N, model)
% given parameters, sample x0
sams = (0.4 + 0.6 * normcdf(thetas(1,:))) .* randn(model.dimension, N);
end