function stats = generate_samples_nlog(model,t, NN)
%{
Generate samples provide model.N samples of {\theta, X_1:t} together with
their ftt probability, theretical probability and ESS
Input: model: structure containing all information till time T
           t: the time step at which all stats are computed
output: stats: structure containing all stats needed
%}
tic
stats.samples = zeros(model.dimension, t+1, model.N*NN);
thetas = zeros(model.td, model.N*NN);
stats.pdf_e = zeros(1, model.N*NN);
stats.pdf_t = zeros(1, model.N*NN);

ftt = model.sirt{t};
[r,f] = eval_irt(ftt, rand(model.d, model.N*NN)); 

switch model.trans
    case 'powertail'
        temp2 = model.L_all(:, :, t) * powinvlag(r, model.pow) + model.mu_all(:, t); % recover samples using L and mu
    case 'truncated'
        temp2 = model.L_all(:, :, t) * r + model.mu_all(:, t); % recover samples using L and mu
end
stats.samples(:, t+1, :) = temp2(model.d - model.dimension*2 + 1 : model.d - model.dimension, :); % generate smaples
stats.samples(:, t, :) = temp2(model.d - model.dimension + 1:model.d, :); % generate smaples

if model.td > 0
%         model.thetas = recovertheta(temp2(1:model.td,:));
    thetas = temp2(1:model.td, :); 
end

switch model.trans
    case 'powertail'
        stats.pdf_e = log(f) + log(powpdflag(powinvlag(r, model.pow), model.pow));
        for k = t-2:(-1):0  % k is time, stored in k+1(second indice) in samples
            ftt = model.sirt{k+1};
            xt = stats.samples(:, k+2, :);
            if model.td > 0
                z = powcdflag(model.L_all(1:model.d-model.dimension,1:model.d-model.dimension, k+1)\([thetas; reshape(xt, model.dimension, model.N*NN)] - model.mu_all(1:model.d-model.dimension,k+1)), model.pow);
            else
                z = powcdflag(model.L_all(1:model.d-model.dimension,1:model.d-model.dimension, k+1)\(reshape(xt, model.dimension, model.N*NN) - model.mu_all(1:model.d-model.dimension,k+1)), model.pow);
            end
            [r,f] = eval_cirt(ftt, z, rand(model.dimension, model.N*NN));
            temp = powinvlag([z ; r], model.pow);
            temp2 = model.L_all(:,:, k+1) * temp + model.mu_all(:, k+1);
            stats.samples(:, k+1, :) = temp2(model.d - model.dimension + 1:model.d, :);
            stats.pdf_e = log(f) + stats.pdf_e - log(powpdflag(temp,model.pow)) - ...
                log(powpdflag(temp(1:model.d - model.dimension, :),model.pow));  
        end
    case 'truncated'
        stats.pdf_e = log(f);
        for k = t-2:(-1):0  % k is time, stored in k+1(second indice) in samples
            ftt = model.sirt{k+1};
            xt = stats.samples(:, k+2, :);
            if model.td > 0
                z = model.L_all(1:model.d-model.dimension,1:model.d-model.dimension, k+1)\([thetas; reshape(xt, model.dimension,[])] - model.mu_all(1:model.d-model.dimension,k+1));
            else
                z = model.L_all(1:model.d-model.dimension,1:model.d-model.dimension, k+1)\(reshape(xt, model.dimension, []) - model.mu_all(1:model.d-model.dimension,k+1));
            end
            [r,f] = eval_cirt(ftt, z, rand(model.dimension, model.N*NN));
            temp = [z ; r];
            temp2 = model.L_all(:,:, k+1) * temp + model.mu_all(:, k+1);
            stats.samples(:, k+1, :) = temp2(model.d - model.dimension + 1:model.d, :);
            stats.pdf_e = log(f) + stats.pdf_e;
        end
end

% calculate ESS
% theoretical pdf
for k = 1:model.N*NN
    dat = stats.samples(:, :, k);
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
w_all = stats.pdf_t - stats.pdf_e;
for k = 1:NN
    w = w_all((k-1)*model.N+1: k*model.N);
    w = exp(w - max(w));
    w = w/nansum(w);
    stats.ess(k) = nansum(w)^2/nansum(w.^2);
end
stats.time_sample = toc;
end

