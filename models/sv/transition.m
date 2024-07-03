function pdfval = transition(~, thetax, ~)
% nx = size(thetax,2);


gamma = normcdf(thetax(1, :))* 0.8 + 0.1;
% ita = sqrt(1-gamma.^2);

X_mu = gamma.* thetax(4, :);
pdfval = normpdf(thetax(3, :) - X_mu);
end