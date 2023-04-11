function stats = generate_samples_dirtlog(model,t, NN)
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
if strcmp(model.coll, 'time1')
    Tthbs = @(x) ftt2Tth(x, model, 1);
    Tthinvbs = @(x) ftt2Tthinv(x, model, 1);
    pdfthbs =  @(x) eval_pdf_check(model.sirt{1}, model.L_all(1: model.td, 1: model.td, 1) \ (x - model.mu_all(1:model.td, 1)), model.poly);
end
L = model.poly.domain(1);
U = model.poly.domain(2);
if t>= 3
    % for time t
    sirt = model.sirt{t};
    [r, f] = eval_irt(sirt, rand(model.d, model.N*NN)); 
    f = log(f) - max(log(f));
    temp = model.L_all(:, :, t) * r + model.mu_all(:, t);
    % recover theta
    thetas = ftt2Tth(temp(1:model.td, :), model, t-1);
    if strcmp(model.coll, 'time1')
        thetas = Tthbs(thetas);
    elseif strcmp(model.coll, 'adapted')
        collind = model.collind(t-1);
        Tthbs = @(u) model.L_coll(1:model.td, 1:model.td, collind) * eval_irt(model.collirt{collind}, model.R(u)) + model.mu_coll(1:model.td, collind);
        thetas = Tthbs(thetas);
    end
    stats.theta_samples = thetas;
    % recover x_k and x_k-1
    X_samples(:, t+1, :) = temp(model.d - model.dimension*2 + 1 : model.d - model.dimension, :);
    X_samples(:, t, :) = ftt2tx([temp(1:model.td, :);temp(model.d - model.dimension + 1:model.d, :)], model, t-1);
    % recover pdf
    
    if strcmp(model.coll, 'NA')
        eval_input_temp = [ftt2Tthinv(thetas, model, t-2); reshape(X_samples(:, t, :), model.dimension, [])] ;
        eval_input = model.L_all(1:model.td + model.dimension,1:model.td + model.dimension, t-1)\...
            (eval_input_temp - model.mu_all(1:model.td + model.dimension, t-1)); % transformation needed on theta to match t-1 input
        if t == 3
            stats.pdf_e = f + log(eval_pdf_check(model.sirt{t-1}, eval_input, model.poly)) - log(model.fu(temp([1:model.td, model.td+model.dimension+1:model.d], :))) +...
            log(eval_pdf_check(model.sirt{t-2}, model.L_all(1:model.td,1:model.td,t-2)\(thetas-model.mu_all(1:model.td,t-2)), model.poly))-log(model.fu(eval_input_temp(1:model.td, :)));
        else
            stats.pdf_e = f + log(eval_pdf_check(model.sirt{t-1}, eval_input, model.poly)) - log(model.fu(temp([1:model.td, model.td+model.dimension+1:model.d], :))) +...
            log(eval_pdf_check(model.reirt{t-2}, model.L_theta(1:model.td,1:model.td,t-2)\(thetas-model.mu_theta(1:model.td,t-2)), model.poly)) - log(model.fu(eval_input_temp(1:model.td, :)));
        end
    elseif strcmp(model.coll, 'time1')
         if t == 3
            eval_input_temp = [ftt2Tthinv(thetas, model, t-2); reshape(X_samples(:, t, :), model.dimension, [])];
            eval_input = model.L_all(1:model.td + model.dimension,1:model.td + model.dimension, t-1)\...
            (eval_input_temp - model.mu_all(1:model.td + model.dimension, t-1)); % transformation needed on theta to match t-1 input
            
            stats.pdf_e = f+log(eval_pdf_check(model.sirt{t-1}, eval_input, model.poly)) - log(model.fu(temp([1:model.td, model.td+model.dimension+1:model.d], :))) +...
            log(eval_pdf_check(model.sirt{t-2}, model.L_all(1:model.td,1:model.td,t-2)\(thetas-model.mu_all(1:model.td,t-2)), model.poly)) - log(model.fu(eval_input_temp(1:model.td, :)));
         else
            eval_input_temp = [ftt2Tthinv(Tthinvbs(thetas), model, t-2); reshape(X_samples(:, t, :), model.dimension, [])];
            eval_input = model.L_all(1:model.td + model.dimension,1:model.td + model.dimension, t-1)\...
            (eval_input_temp - model.mu_all(1:model.td + model.dimension, t-1)); % transformation needed on theta to match t-1 input
           
            stats.pdf_e = f + log(eval_pdf_check(model.sirt{t-1}, eval_input, model.poly))- log(model.fu(temp([1:model.td, model.td+model.dimension+1:model.d], :))) +...
            log(eval_pdf_check(model.reirt{t-2}, model.L_theta(1:model.td,1:model.td,t-2)\(Tthinvbs(thetas)-model.mu_theta(1:model.td,t-2)), model.polycoll)) -...
            log(model.fu(eval_input_temp(1:model.td, :))) + log(pdfthbs(thetas)) - log(model.fu(Tthinvbs(thetas)));
         end
    elseif strcmp(model.coll, 'adapted') 
        if t == 3
            eval_input_temp = [ftt2Tthinv(thetas, model, t-2); reshape(X_samples(:, t, :), model.dimension, [])];
            eval_input = model.L_all(1:model.td + model.dimension,1:model.td + model.dimension, t-1)\...
            (eval_input_temp - model.mu_all(1:model.td + model.dimension, t-1)); % transformation needed on theta to match t-1 input
            
            stats.pdf_e = f + log(eval_pdf_check(model.sirt{t-1}, eval_input, model.poly)) - log(model.fu(temp([1:model.td, model.td+model.dimension+1:model.d], :))) +...
            log(eval_pdf_check(model.sirt{t-2}, model.L_all(1:model.td,1:model.td,t-2)\(thetas-model.mu_all(1:model.td,t-2)), model.poly)) - log(model.fu(eval_input_temp(1:model.td, :)));
        else
            collind = model.collind(t-2);
            Tthinvbs = @(x) model.Rinv(eval_rt_check(model.collirt{collind}, model.L_coll(1:model.td,1:model.td, collind)\(x-model.mu_coll(1:model.td, collind)), model.polycoll));
            pdfthbs = @(x) eval_pdf_check(model.collirt{collind},model.L_coll(1:model.td,1:model.td, collind) \ (x - model.mu_coll(1:model.td, collind)), model.polycoll);
            
            eval_input_temp = [ftt2Tthinv(Tthinvbs(thetas), model, t-2); reshape(X_samples(:, t, :), model.dimension, [])];
            eval_input = model.L_all(1:model.td + model.dimension,1:model.td + model.dimension, t-1)\...
            (eval_input_temp - model.mu_all(1:model.td + model.dimension, t-1)); % transformation needed on theta to match t-1 input
           
            stats.pdf_e = f + log(eval_pdf_check(model.sirt{t-1}, eval_input, model.poly)) - log(model.fu(temp([1:model.td, model.td+model.dimension+1:model.d], :))) +...
            log(eval_pdf_check(model.reirt{t-2}, model.L_theta(1:model.td,1:model.td,t-2)\(Tthinvbs(thetas)-model.mu_theta(1:model.td,t-2)), model.polycoll)) -...
            log(model.fu(eval_input_temp(1:model.td, :))) + log(pdfthbs(thetas)) -log(model.fu(Tthinvbs(thetas)));
         end
    else 
        error('no collapse precondition option')
    end


    % from X_t-2 to 2

    for k = t-2:-1:2 % k is time, stored in k+1(second indice) in samples
        sirt = model.sirt{k+1};
        switch model.coll
            case 'NA'
                cond = model.L_all(1:model.d-model.dimension,1:model.d-model.dimension, k+1) \...
            ([ftt2Tthinv(thetas, model, k) ;reshape(X_samples(:, k+2, :), model.dimension, [])] - model.mu_all(1:model.d-model.dimension, k+1));
            case 'time1'
                cond = model.L_all(1:model.d-model.dimension,1:model.d-model.dimension, k+1) \...
            ([ftt2Tthinv(Tthinvbs(thetas), model, k) ;reshape(X_samples(:, k+2, :), model.dimension, [])] - model.mu_all(1:model.d-model.dimension, k+1));
            case 'adapted'
                collind = model.collind(k);
                Tthinvbs = @(x) model.Rinv(eval_rt(model.collirt{collind}, model.L_coll(1:model.td,1:model.td, collind)\(x-model.mu_coll(1:model.td, collind))));
                cond = model.L_all(1:model.d-model.dimension,1:model.d-model.dimension, k+1) \...
            ([ftt2Tthinv(Tthinvbs(thetas), model, k) ;reshape(X_samples(:, k+2, :), model.dimension, [])] - model.mu_all(1:model.d-model.dimension, k+1));
        end
        ind = all((cond <= U & cond >= L), 1); % flag to check the range
        dex = find(indall == 1);
        indall(dex(ind == 0)) = 0;
%         fprintf('%d samples are missing at step %d\n', length(rg_flag)-nnz(rg_flag), k)
        cond = cond(:, ind); % points out of range are omitted so that no NaN in future computations
        
        [r,f] = eval_cirt(sirt, cond, rand(model.dimension, size(cond, 2)));
        f = log(f) - max(log(f));
        temp = model.L_all(:, :, k+1) * [cond; r] + model.mu_all(:, k+1);
        X_samples = X_samples(:, :, ind);
        X_samples(:, k+1, :) = ftt2tx([temp(1:model.td, :);temp(model.d - model.dimension + 1:model.d, :)], model, k);
        thetas = thetas(:, ind);
        if strcmp(model.coll, 'NA')
            % transformation needed on theta to match t-1 input
            eval_input = model.L_all(1:model.td + model.dimension,1:model.td + model.dimension, k)\...
            ([ftt2Tthinv(thetas, model, k-1); reshape(X_samples(:, k+1, :), model.dimension, [])] - model.mu_all(1:model.td + model.dimension, k));
        elseif strcmp(model.coll, 'time1')
            if k == 2
            eval_input = model.L_all(1:model.td + model.dimension,1:model.td + model.dimension, k)\...
            ([ftt2Tthinv(thetas, model, k-1); reshape(X_samples(:, k+1, :), model.dimension, [])] - model.mu_all(1:model.td + model.dimension, k));
            else
            eval_input = model.L_all(1:model.td + model.dimension,1:model.td + model.dimension, k)\...
            ([ftt2Tthinv(Tthinvbs(thetas), model, k-1); reshape(X_samples(:, k+1, :), model.dimension, [])] - model.mu_all(1:model.td + model.dimension, k));
            end
        elseif strcmp(model.coll, 'adapted')
            if k == 2
                Tthinvbs = @(x) x;
            else
                collind = model.collind(k-1);
                Tthinvbs = @(x) model.Rinv(eval_rt(model.collirt{collind}, model.L_coll(1:model.td,1:model.td, collind)\(x-model.mu_coll(1:model.td, collind))));
            end
            eval_input = model.L_all(1:model.td + model.dimension,1:model.td + model.dimension, k)\...
            ([ftt2Tthinv(Tthinvbs(thetas), model, k-1); reshape(X_samples(:, k+1, :), model.dimension, [])] - model.mu_all(1:model.td + model.dimension, k));
        else
            error('no collapse precondition option')
        end
        
        stats.pdf_e = stats.pdf_e(ind) + f + log(eval_pdf_check(model.sirt{k}, eval_input, model.poly))...
            - log(model.fu(temp([1:model.td, model.td+model.dimension+1:model.d], :)))...
            + log(model.fu(temp(1:model.td, :)))- log(eval_pdf_check(model.sirt{k}, eval_input(1:model.td, :), model.poly));
%         stats.pdf_e = stats.pdf_e + f + log(eval_pdf_check(model.sirt{k}, eval_input, model.poly))...
%             - log(model.fu(temp([1:model.td, model.td+model.dimension+1:model.d], :)))...
%             + log(model.fu(temp(1:model.td, :)))- log(eval_pdf_check(model.sirt{k}, eval_input(1:model.td, :), model.poly));
    end
    % sample for X_1
    sirt = model.sirt{2};
    cond = model.L_all(1:model.d-model.dimension,1:model.d-model.dimension, 2) \...
        ([ftt2Tthinv(thetas, model, 1) ;reshape(X_samples(:, 3, :), model.dimension, [])] - model.mu_all(1:model.d-model.dimension, 2));
    ind = all((cond <= U & cond >= L), 1); % flag to check the range
    dex = find(indall == 1);
    indall(dex(ind == 0)) = 0;
    cond = cond(:, ind);
    X_samples = X_samples(:, :, ind);
    thetas = thetas(:, ind);
    [r,f] = eval_cirt(sirt, cond, rand(model.dimension, size(cond, 2)));
    f = log(f) - max(log(f));
    temp = model.L_all(:, :, 2) * [cond; r] + model.mu_all(:, 2);
    X_samples(:, 2, :) = ftt2tx([temp(1:model.td, :);temp(model.d - model.dimension + 1:model.d, :)], model, 1);
    eval_input = model.L_all(1:model.td + model.dimension,1:model.td + model.dimension, 1)\...
    ([thetas; reshape(X_samples(:, 2, :), model.dimension, [])] - model.mu_all(1:model.td + model.dimension, 1)); % transformation needed on theta to match t-1 input
%     stats.pdf_e = stats.pdf_e + f + log(eval_pdf_check(model.sirt{1}, eval_input, model.poly))...
%         - log(model.fu(temp([1:model.td, model.td+model.dimension+1:model.d], :)))...
%         + log(model.fu(temp(1:model.td, :))) - log(eval_pdf_check(model.sirt{1}, eval_input(1:model.td, :), model.poly));
    stats.pdf_e = stats.pdf_e(ind) + f + log(eval_pdf_check(model.sirt{1}, eval_input, model.poly))...
        - log(model.fu(temp([1:model.td, model.td+model.dimension+1:model.d], :)))...
        + log(model.fu(temp(1:model.td, :))) - log(eval_pdf_check(model.sirt{1}, eval_input(1:model.td, :), model.poly));


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
elseif t == 2
    sirt = model.sirt{2};
    [r,f] = eval_irt(sirt, rand(model.d, model.N*NN));
    f = log(f) - max(log(f));
    temp = model.L_all(:, :, 2) * r + model.mu_all(:, 2);
    X_samples(:, 2, :) = ftt2tx([temp(1:model.td, :);temp(model.d - model.dimension + 1:model.d, :)], model, 1);
    X_samples(:, 3, :) = temp(model.d-2*model.dimension+1:model.d-model.dimension, :);
    thetas = ftt2Tth(temp(1:model.td, :), model, 1);

    eval_input = model.L_all(1:model.td + model.dimension,1:model.td + model.dimension, 1)\...
    ([thetas; reshape(X_samples(:, 2, :), model.dimension, [])] - model.mu_all(1:model.td + model.dimension, 1));
    stats.pdf_e = f + log(eval_pdf_check(model.sirt{1}, eval_input, model.poly))...
        - log(model.fu(temp([1:model.td, model.td+model.dimension+1:model.d], :)));

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

