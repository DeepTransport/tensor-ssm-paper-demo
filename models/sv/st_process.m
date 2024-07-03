function X_new = st_process(~, thetax, ~)
nx = size(thetax, 2);



gamma = normcdf(thetax(1, :))*0.8+0.1;
% ita = sqrt(1-gamma.^2);

X_new = gamma.* thetax(3, :) + 1 .* randn(1, nx);
end