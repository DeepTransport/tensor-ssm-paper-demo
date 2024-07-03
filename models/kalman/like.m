function pdf = like(model, thetax, t)
y = model.Y(:, t);
a = .4 + .6 * normcdf(thetax(1, :));
d = .4 + .6 * normcdf(thetax(2, :));
% n = size(thetax,2);

yd = y./d;
mu = model.pre.C * (d.^(-1).* thetax(model.d+1:model.d+model.m, :));

if size(y, 1) == 1
    pdf = normpdf(yd, mu)./(d.^model.n);
else
    pdf = mvnpdf(yd', mu')'./(d.^model.n);
end
% pdf = zeros(1,n);
% for k = 1:n
%     d = normcdf(thetax(2,k)) + 0.5;
%     pdf(k) = mvnpdf(y, model.pre.C * thetax(model.d+1:model.d + model.m,k), d^2*eye(model.n));
% end
end