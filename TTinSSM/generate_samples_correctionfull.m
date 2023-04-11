function stats = generate_samples_correctionfull(model,t, NN)
%{
Generate samples provide model.N samples of {\theta, X_1:t} together with
their ftt probability, theretical probability and ESS
Input: model: structure containing all information till time T
           t: the time step at which all stats are computed
output: stats: structure containing all stats needed
%}
tic
indall = true(1, model.N*NN);
X_samples = zeros(model.dimension, t+1, model.N*NN);

L = model.poly.domain(1);
U = model.poly.domain(2);
if t>= 2
    % for time t
    sirt = model.sirt{t};
    [r, f] = eval_irt(sirt, rand(model.d, model.N*NN)); 
    f = log(f) - max(log(f));
    temp = r;
    % recover theta
    thetax = model.L_theta(:,:, t-1) * eval_irt(model.reirt{t-1}, ...
        model.R(temp)) + model.mu_theta(:,t-1);
    thetas = thetax(1:model.td, :);
    stats.theta_samples = thetas;
    % recover x_k and x_k-1
    X_samples(:, t+1, :) = thetax(model.d - model.dimension*2 + 1 : model.d - model.dimension, :);
    X_samples(:, t, :) = thetax(model.td + model.dimension +1:model.d, :);
    % recover pdf
    
    if strcmp(model.coll, 'NA')
        L_invtx = model.L_theta(:,:,t-1)\([thetas; reshape(X_samples(:, t+1, :), model.dimension, []); reshape(X_samples(:, t, :), model.dimension, [])]-model.mu_theta(:, t-1));
        T_invtx = model.Rinv(eval_rt(model.reirt{t-1}, L_invtx));  % T^{-1} (theta, xT-1)
        stats.pdf_e = f + log(eval_pdf_check(model.reirt{t-1}, L_invtx , model.poly))- log(model.fu(T_invtx));
    else 
        error('no collapse precondition option')
    end


    % from X_t-2 to 2

    for k = t-2:-1:1 % k is time, stored in k+1(second indice) in samples
        sirt = model.sirt{k+1};
        switch model.coll
            case 'NA'
                T_invtxp = model.Rinv(eval_rt(model.reirt{k}, model.L_theta(1:model.td+model.dimension, 1:model.td+model.dimension, k)\...
                    ([thetas; reshape(X_samples(:, k+2, :), model.dimension, [])] - model.mu_theta(1:model.td+model.dimension, k))));
                cond = T_invtxp;
        end
        ind = all((cond <= U & cond >= L), 1); % flag to check the range
        dex = find(indall == 1);
        indall(dex(ind == 0)) = 0;
%         fprintf('%d samples are missing at step %d\n', length(rg_flag)-nnz(rg_flag), k)
        T_invtxp = T_invtxp(:, ind);
        cond = cond(:, ind); % points out of range are omitted so that no NaN in future computations
        
        [r,f] = eval_cirt(sirt, cond, rand(model.dimension, size(cond, 2)));
        f = log(f) - max(log(f));
        temp = [cond; r];
        temp = model.L_theta(:,:, k) * eval_irt(model.reirt{k}, ...
            model.R(temp)) + model.mu_theta(:,k);
        X_samples = X_samples(:, :, ind);
        X_samples(:, k+1, :) = temp(model.td+model.dimension+1:model.d, :);
        thetas = thetas(:, ind);
        if strcmp(model.coll, 'NA')
            L_invtx = model.L_theta(:,:,k)\([thetas;reshape(X_samples(:, k+2, :), model.dimension, []); reshape(X_samples(:, k+1, :), model.dimension, [])]-model.mu_theta(:, k));
            T_invtx = model.Rinv(eval_rt(model.reirt{k}, L_invtx));  % T^{-1} (theta, xT-1)
            stats.pdf_e = stats.pdf_e(ind) + f + log(eval_pdf_check(model.reirt{k}, L_invtx, model.polycoll))...
            - log(eval_pdf_check(model.reirt{k}, L_invtx(1:model.td+model.dimension, :), model.polycoll))...
            + log(model.fu(T_invtx(1:model.td+model.dimension,:))) - log(model.fu(T_invtx));
        else
            error('no collapse precondition option')
        end
    end
  

    % sample for X_0
    sirt = model.sirt{1};
    cond = model.L_all(1:model.d-model.dimension,1:model.d-model.dimension, 1) \...
            ([thetas ;reshape(X_samples(:, 2, :), model.dimension, [])] - model.mu_all(1:model.d-model.dimension, 1));
    ind = all((cond <= U & cond >= L), 1); % flag to check the range
    dex = find(indall == 1);
    indall(dex(ind == 0)) = 0;
    cond = cond(:, ind);
    X_samples = X_samples(:, :, ind);
    thetas = thetas(:, ind);
    [r,f] = eval_cirt(sirt, cond, rand(model.dimension, size(cond, 2)));
    f = log(f) - max(log(f));
    temp = model.L_all(:, :, 1) * [cond; r] + model.mu_all(:, 1);
    X_samples(:, 1, :) = temp(model.d-model.dimension+1: model.d, :);
%     stats.pdf_e = stats.pdf_e + f;
    stats.pdf_e = stats.pdf_e(ind) + f;


    stats.theta_samples = thetas;
    stats.X_samples = X_samples; 
    
elseif t == 1
    sirt = model.sirt{1};
    [r,f] = eval_irt(sirt, rand(model.d, model.N*NN));
    f = log(f) - max(log(f));
    temp = model.L_all(:, :, 1) * r + model.mu_all(:, 1);
    X_samples(:, 2, :) = temp(model.d-2*model.dimension+1:model.d-model.dimension, :);
    X_samples(:, 1, :) = temp(model.d-model.dimension+1: model.d, :);
    stats.pdf_e = f;
    thetas = temp(1:model.td, :);
    stats.theta_samples = thetas;
    stats.X_samples = X_samples;   
else
    error('wrong time input')
end
    
    
    
% calculate ESS
% theoretical pdf
stats.pdf_t = zeros(1, length(stats.pdf_e));
for k = 1:length(stats.pdf_e)
    dat = X_samples(:, :, k);
    X0 = [repmat(thetas(:, k), 1, size(dat,2)-1); dat(:, 1:end-1)];

    X1 = X0;
    X1(model.td + 1 : model.td+model.dimension, :) = dat(:, 2:end);
    X1(model.td+model.dimension+1 : model.d, :) = X0(model.td+1 : model.td+model.dimension, :);
    stats.pdf_t(k) = log(prior(X0(:, 1), model)); % this is the prior of X_0 and theta
    for l = 1:size(dat,2)-1
        stats.pdf_t(k) = stats.pdf_t(k) + log(transition(X1(: ,l), model)) + log(like(X1(:, l), model, model.Y(:, l)));
    end
end


stats.ess = zeros(NN, 1);
w_all = zeros(1, model.N*NN)-inf;
w_all(indall) = stats.pdf_t - stats.pdf_e;
w_all(~indall) = -inf;
w_all(isinf(w_all)) = -inf;
w_all = w_all - max(w_all(indall));

for k = 1:NN
    w = w_all((k-1)*model.N+1: k*model.N);
    w = exp(w - max(w(~isinf(w))));
    w = w/nansum(w);
    stats.ess(k) = nansum(w)^2/nansum(w.^2);
end
stats.time_sample = toc;
end

