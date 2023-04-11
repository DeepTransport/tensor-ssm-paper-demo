function value = ftt2tx(x, model, t)
% the function is part of the HMM_FTT project
% the function creates the function handle transport map for x at time t
sirt = model.sirt{t};
value = model.xmat * (model.L_all(1: model.td+model.dimension, 1: model.td+model.dimension, t)...
    * eval_irt(sirt, model.R(x)) + model.mu_all(1: model.td+model.dimension, t));
end

