function [Tth, Tthinv] = ftt2th(model, t)
% the function is part of the HMM_FTT project
% the function creates the function handle transport map for theta at time t
if t == 1
    sirt = model.sirt{1};
    Tth = @(u) model.L_all(1: model.td, 1: model.td, t) * eval_irt(sirt, model.R(u)) + model.mu_all(1: model.td, t);
    Tthinv = @(x) model.Rinv(eval_rt(sirt, model.L_all(1: model.td, 1: model.td, t)\(x-model.mu_all(1: model.td, t)))); 
else
    sirt = model.reirt{t}; 
    Tth = @(u) model.L_theta(:,:,t) * eval_irt(sirt, model.R(u)) + model.mu_theta(:,t);
    Tthinv = @(x) model.Rinv(eval_rt(sirt, model.L_theta(:,:,t)\(x-model.mu_theta(:,t))));
end
end

