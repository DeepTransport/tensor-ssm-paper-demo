function pdf = priorpdf(model, thetax)
a = .4 + .6*normcdf(thetax(1, :));
n = size(thetax,2);
pdf = zeros(1,n);
for k = 1:n
    pdf(k) = mvnpdf(thetax(1:model.d, k)', zeros(model.d, 1)', eye(model.d));
    pdf(k) = pdf(k) * mvnpdf(thetax(1+model.d:model.d+model.m, k)', zeros(model.m, 1)',  eye(model.m));
% pdf(k) = pdf(k) * mvnpdf(thetax(1+model.d:end, k)', zeros(model.m, 1)',eye(model.m));
end
end