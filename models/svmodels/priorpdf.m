function pdfval = priorpdf(model, thetax)

[gamma, tau, mu, phi, a, delta, nu1, nu2] = ftt2true(thetax(1:model.d, :), model);


invgamapdf = (1./(tau.^2)).^2 .* exp(- 0.005 ./ (tau.^2)) * 0.005;

pdfx = normpdf(thetax(model.d+1, :), mu./tau, 1./sqrt(1-gamma.^2));
pdf12 = betapdf((normcdf(thetax(1,:))+1)/2, 20, 1.5).* normpdf(thetax(1,:)) .* ...
    invgamapdf.*tau.^2;

if isnan(model.pre.para(3))
    pdf3 = normpdf(thetax(3, :), 0, sqrt(5));
    ind = 1;
else
    pdf3 = ones(1, size(thetax, 2));
    ind = 0;
end

if isnan(model.pre.para(4)) && isnan(model.pre.para(5))
    pdf45 = prod(normpdf(thetax(3+ind:4:ind, :), 0, sqrt(2)));
    ind = ind + 2;
elseif isnan(model.pre.para(4)) || isnan(model.pre.para(5))
    pdf45 = normpdf(thetax(3+ind, :), 0, sqrt(2));
    ind = ind + 1;
else
    pdf45 = ones(1, size(thetax, 2));
    ind = ind + 0;
end

pdf678 = prod(normpdf(thetax(3+ind:model.d, :), 0, 1), 1);


pdfval = pdf12.*pdf3.*pdf45.*pdf678.*pdfx;
pdfval(isnan(pdfval)) = 0;

% pdfval = normpdf(thetax(1,:)).*normpdf(thetax(2,:)).*normpdf(thetax(3,:)).*normpdf(thetax(4,:)).*normpdf(thetax(5,:), thetax(2,:));
end