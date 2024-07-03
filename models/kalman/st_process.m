function X_new = st_process(model, thetax, ~)

n = size(thetax,2);
a = .4 + .6*normcdf(thetax(1, :));
X_new = sqrt(1-a.^2).*  thetax(model.d+1 : model.d + model.m, :) + a .* randn(model.m, n)  ;

end