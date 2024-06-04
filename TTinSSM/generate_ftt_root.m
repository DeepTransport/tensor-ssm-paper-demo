function model = generate_ftt_root(model)

% sam_old = priorsam(model);
fun = @(X) prior(X, model);

switch model.trans       
    case 'truncated'
        for k = 1 : model.steps
            fprintf('this is ftt build up at time %d\n', k)
            tic        
            L = eye(model.td + 2*model.dimension);
            mu = zeros(model.td + 2*model.dimension, 1);
            L1 = L(1:model.td + model.dimension, 1:model.td + model.dimension);
            model.L_all(:, :, k) = L;
            mu1 = mu(1:model.td + model.dimension);
            model.mu_all(:, k) = mu;
            
            fun = @(z) fun(model.cmat*(L*z + mu)).*transition(L*z + mu, model).*...
                 like(L*z + mu, model, model.Y(:,k)); % update
            if k == 1
                ftt = FTT(fun, model.d, model.poly, model.opt);
            else 
                ftt = FTT(fun, model.d, ftt);
            end

             
            fun = @(x) eval(int_block(ftt, model.td + model.dimension+1:model.d), L1\(x-mu1)); % input is inf       
            model.time_filter(k) = toc;  
        end                        
end