function Y = ob_process(model, thetax, ~)
% this is only used to generate observations so use true theta
% X only contains states
th = normcdf(model.theta) * 0.6 + 0.4;
Y = ( model.pre.C * thetax(model.d+1:end,:) )+ th(2) * randn(model.n, size(thetax(model.d+1:end,:), 2));
end