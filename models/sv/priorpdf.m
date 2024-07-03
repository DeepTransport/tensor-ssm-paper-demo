function pdfval = priorpdf(model, thetax)

gamma = normcdf(thetax(1, :)) *0.8+0.1;
% mu = thetax(2, :);
% tau2 = normcdf(thetax(2, :));
% phi = thetax(4, :);

pdfval =  normpdf(thetax(model.d+1, :), 0, 1./sqrt(1-gamma.^2)) .* ...
    normpdf(thetax(1, :)) .*...
    normpdf(thetax(2, :));
pdfval(isnan(pdfval)) = 0;
end