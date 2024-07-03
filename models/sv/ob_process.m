function Y = ob_process(model, thetax, ~)
% this is only used to generate observations so use true theta
% X only contains states
% rng(1)
Y = (normcdf(model.theta(2))*0.8+0.1) * exp(0.5 * thetax(model.d+1:end, :)) * randn(1,size(thetax, 2));
end