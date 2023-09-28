function model = generate_ftt_nl(model)

sam_old = priorsam(model);
fun = @(X) prior(X, model);

switch model.trans
    case 'powertail'
        for k = 1 : model.steps
            fprintf('this is ftt build up at time %d\n', k)
            tic
            sam_temp = st_process(sam_old, model);
            sam_new = sam_old;
            sam_new(model.td + 1 : model.td + model.dimension, :) = sam_temp;
            sam_new(model.td + model.dimension+1 : model.d, :) = sam_old(model.td + 1:model.td + model.dimension, :);
            w = like(sam_new, model, model.Y(:, k));
            sam_new = datasample(sam_new', model.sam_size,'Weights',w);
            sam_new = sam_new';
            stepess = sum(w)^2/sum(w.^2)/model.sam_size*100;
            model.stepess(k) = stepess;
            fprintf('ESS at step %d is %.2f \n', k, stepess)
            
            mu = mean(sam_new, 2);
            L = chol(cov(sam_new') + 1e-8 * eye(size(sam_new,1)))'; % linear transformation
            [~, L_quantile] = trans(L\(sam_new-mu), model.quantile, model.pow); % this step is to scale transformation
            L = L*L_quantile;
            % get block matrices
            L1 = L(1:model.td + model.dimension, 1:model.td + model.dimension);
            model.L_all(:, :, k) = L;
            mu1 = mu(1:model.td + model.dimension);
            model.mu_all(:, k) = mu;
            
            fun = @(z) nan2zero(fun(model.cmat*(L*powinvlag(z, model.pow) + mu)).*transition(L*powinvlag(z, model.pow) + mu, model).*...
                 like(L*powinvlag(z, model.pow) + mu, model, model.Y(:,k))./powpdflag(powinvlag(z, model.pow), model.pow)); % update
            
            if k == 1
                ftt = SIRT(fun, model.d, model.poly, model.opt, 'sample_x', ...
                    powcdflag(L\(sam_new(:, 1:model.sam_size/2)-mu), model.pow), 'debug_x', powcdflag(L\(sam_new(:, model.sam_size/2:model.sam_size)-mu), model.pow));
            else
                ftt = SIRT(fun, model.d, ftt, 'sample_x', ...
                    powcdflag(L\(sam_new(:, 1:model.sam_size/2)-mu), model.pow), 'debug_x', powcdflag(L\(sam_new(:, model.sam_size/2:model.sam_size)-mu), model.pow));
            end
            model.sirt{k} = ftt;      
            sam_old = eval_irt(ftt, rand(model.d - model.dimension, model.sam_size));
            sam_old = L1 * powinvlag(sam_old,model.pow) + mu1; % generate samples     
            fun = @(x) eval_pdf(ftt, powcdflag(L1\(x-mu1),model.pow)).*powpdflag(L1\(x-mu1),model.pow); % input is inf     
            model.time_filter(k) = toc;  
        end
        
    case 'truncated'
        for k = 1 : model.steps
            fprintf('this is ftt build up at time %d\n', k)
            tic
            sam_temp = st_process(sam_old, model);
            sam_new = sam_old;
            sam_new(model.td+1:model.td+model.dimension, :) = sam_temp;
            sam_new(model.td+model.dimension+1:model.d, :) = sam_old(model.td+1:model.td+model.dimension, :);
            w = like(sam_new, model, model.Y(:, k));
            sam_new = datasample(sam_new', model.sam_size,'Weights',w);
            sam_new = sam_new';
            stepess = sum(w)^2/sum(w.^2)/model.sam_size*100;
            model.stepess(k) = stepess;
            fprintf('ESS at step %d is %.2f \n', k, stepess)
           
             [L, mu] = ab(sam_new, model.sdn);  % this step is to scale transformation
%             mu = mean(sam_new, 2);
%             L = chol(cov(sam_new') + 1e-8 * eye(model.d))';
            
            % get block matrices
            L1 = L(1:model.td + model.dimension, 1:model.td + model.dimension);
            model.L_all(:, :, k) = L;
            mu1 = mu(1:model.td + model.dimension);
            model.mu_all(:, k) = mu;
            
            fun = @(z) fun(model.cmat*(L*z + mu)).*transition(L*z + mu, model).*...
                 like(L*z + mu, model, model.Y(:,k)); % update
            if k == 1
                ftt = SIRT(fun, model.d, model.poly, model.opt, 'sample_x',...
                    L\(sam_new(:, 1:model.sam_size/2)-mu), 'debug_x', L\(sam_new(:, model.sam_size/2:model.sam_size)-mu));
            else 
                ftt = SIRT(fun, model.d, ftt, 'sample_x',...
                    L\(sam_new(:, 1:model.sam_size/2)-mu), 'debug_x', L\(sam_new(:, model.sam_size/2:model.sam_size)-mu));
            end
            [sam, psam] = eval_irt(ftt, rand(model.d, 1e4));
            w = fun(sam)./psam;
            model.fullESS(k) = sum(w)^2/sum(w.^2);

            model.sirt{k} = ftt;
            sam_old = eval_irt(ftt, rand(model.d - model.dimension, model.sam_size));
            sam_old = L1 * sam_old + mu1; % generate samples
             
            fun = @(x) eval_pdf(ftt, L1\(x-mu1)); % input is inf       
            model.time_filter(k) = toc;  
        end                        
end