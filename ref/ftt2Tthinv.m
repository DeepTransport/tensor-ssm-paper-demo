function value = ftt2Tthinv(x, model, t)
% the function is part of the HMM_FTT project
% the function creates the function handle transport map for theta at time t
if t == 1
    sirt = model.sirt{1};
    value = model.Rinv(eval_rt(sirt, model.L_all(1: model.td, 1: model.td, t)\(x-model.mu_all(1: model.td, t)))); 
else
    sirt = model.reirt{t}; 
    value = model.Rinv(eval_rt_check(sirt, model.L_theta(:,:,t)\(x-model.mu_theta(:,t)), model.poly));
end
end

