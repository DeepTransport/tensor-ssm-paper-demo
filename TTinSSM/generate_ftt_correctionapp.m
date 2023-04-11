 function model = generate_ftt_correctionapp(model)
%  GENERATE_FTT Given the model, generate ftt in every step.
%  every step including correction ratio. The collapse step is for theta
%  and x_t together
model.mu_theta = zeros(model.d, model.steps);
model.L_theta = zeros(model.d, model.d, model.steps);
% define some transition matrix
% comment order [theta, x1, x2]
rt = sparse([eye(model.td, model.td), zeros(model.td, model.dimension * 2)]); % transition matrix to get theta
r1 = sparse([zeros(model.dimension, model.td), eye(model.dimension),zeros(model.dimension,model.dimension)]); % transition matrix to get x1
r2 = sparse([zeros(model.dimension, model.td), zeros(model.dimension,model.dimension), eye(model.dimension)]);
rt2 = sparse(model.td+model.dimension, model.d);
rt2(1:model.td, 1:model.td) = eye(model.td); % transition matrix to get theta and x2
rt2(model.td+1:model.td+model.dimension, model.d-model.dimension+1:model.d) = eye(model.dimension);  % transition matrix to get theta and x2
s1 = sparse([eye(model.td), zeros(model.td, model.dimension)]); % get theta from [theta, x2]
s2 = sparse([zeros(model.dimension, model.td), eye(model.dimension)]); % get x2 from [theta, x2]


tic
sam_old = priorsam(model);
fun = @(X) prior(X, model);
% at time 1
% regularization begins
sam_temp = st_process(sam_old, model);
sam_new = sam_old;
sam_new(model.td + 1 : model.td + model.dimension, :) = sam_temp;
sam_new(model.td + model.dimension+1 : model.d, :) = sam_old(model.td + 1:model.td + model.dimension, :);
w = like(sam_new, model, model.Y(:, 1));
sam_new = datasample(sam_new', model.sam_size,'Weights',w);
sam_new = sam_new';

L = model.amplify*diag(sqrt(var(sam_new, 0, 2)));
mu = mean(sam_new, 2);
if strcmp(model.name, 'kalman')
    L = chol(cov(sam_new') + 1e-8 * eye(size(sam_new,1)))';
end

% store block matrices
model.L_all(:, :, 1) = L;
model.mu_all(:, 1) = mu;
% regularization ends

fun = @(x) fun(model.cmat*(L*x+mu)).*transition((L*x+mu), model).*...
        like((L*x+mu), model, model.Y(:,1)); % update
sirt = SIRT(fun, model.d, model.poly, model.optinit);
model.sirt{1} = sirt;

% new samples
sam_uni = rand(model.d - model.dimension, model.sam_size);
sam_u = model.Rinv(sam_uni); % this will be the input samples to next iteration(used to calculate L) 
sam_old = L(1: model.td+model.dimension, 1: model.td+model.dimension) * eval_irt(sirt, sam_uni) + mu(1: model.td+model.dimension);

sam_temp = st_process(sam_old, model);
sam_new = sam_old;
sam_new(model.td+1:model.td+model.dimension, :) = sam_temp;
sam_new(model.td+model.dimension+1:model.d, :) = sam_old(model.td+1:model.td+model.dimension, :);
w = like(sam_new, model, model.Y(:, 2));
stepess = sum(w)^2/sum(w.^2)/model.sam_size*100;
model.stepess(1) = stepess;

fprintf('precollapse at step %d is %.2f \n', 1, stepess)
sam_new = datasample(sam_new', model.sam_size,'Weights',w);
sam_new = sam_new';


mu_theta = mean(sam_new, 2);
L_theta = model.amplify*diag(sqrt(var(sam_new, 0, 2)));
if model.linear == 1
    Sigma = (model.c * cov(sam_new')^(-1) + (1- model.c) * eye(size(sam_new,1)))^(-1);
    L_theta = chol(Sigma + 1e-8 * eye(size(sam_new,1)))';
end

model.L_theta(:, :, 1) = L_theta;
model.mu_theta(:, 1) = mu_theta;

Lx = @(x) L_theta*x+mu_theta;
pdffull = @(x) (nan2zero(eval_pdf(sirt, L(1:model.td+model.dimension,1:model.td+model.dimension)\(rt2*Lx(x) - mu(1:model.td+model.dimension)))) .* ...
    transition(Lx(x), model).*like(Lx(x), model, model.Y(:,2))) ./model.fu(x);
% const = pdffull(zeros(model.d, 1))*model.c;
pdffull = @(x) pdffull(x) .^ model.c.*model.fu(x);

pdfth_app = pdffull;

debug_x = L_theta\(sam_new-mu_theta);

fprintf('\n collapse for Tth at time %d\n', 1)

reirt = SIRT(pdffull, model.d, model.polycoll, model.optcoll, 'sample_x', debug_x(:, 1:model.sam_size/2), 'debug_x', debug_x(:, model.sam_size/2:model.sam_size));

model.reirt{1} = reirt;

pdfthold = @(x) nan2zero(eval_pdf_check(sirt, L(1:model.td+model.dimension,1:model.td+model.dimension)\(x-mu(1:model.td+model.dimension)), model.poly));
That = @(u) L_theta * eval_irt(reirt, model.R(u)) + mu_theta;
Thatinv = @(x) model.Rinv(eval_rt(reirt, L_theta\(x-mu_theta)));
Thatinv2 = @(x) model.Rinv(eval_rt(reirt, L_theta(1:model.td+model.dimension, 1:model.td+model.dimension)\(x-mu_theta(1:model.td+model.dimension))));
pdfth = @(x) nan2zero(eval_pdf_check(reirt, L_theta\(x-mu_theta), model.polycoll));
pdfth2 = @(x) nan2zero(eval_pdf_check(reirt, L_theta(1:model.td+model.dimension, 1:model.td+model.dimension)\(x-mu_theta(1:model.td+model.dimension)), model.polycoll));

thetasams = eval_irt(reirt, rand(model.d, 1000));
w = pdfth_app(thetasams)./eval_pdf(reirt, thetasams);
model.collESS(1) = sum(w)^2/sum(w.^2)/10;
sam_L = Thatinv(sam_new);
% build T_{1}^\theta

model.time_filter(1) = toc;


for k = 2 : model.steps
    tic
    % regularization begins
    
    
    sam_L(:, any(isnan(sam_L))) = [];
    
    mu = zeros(model.d,1);
    L = eye(model.d);  
    model.L_all(:,:,k) = L;
    model.mu_all(:,k) = mu;

    fprintf('this is ftt build up at time %d\n', k)
    fun = @(u) funintoapp(u, mu, L, That, rt2, pdfth, pdfthold, model, k);
    sirtold = sirt;
    debug_x = L\(sam_L-mu);
    location = floor(size(debug_x,2)/2);
    if k == 2 
        sirt = SIRT(fun, model.d, model.poly, model.opt, 'sample_x', debug_x(:, 1:location), 'debug_x', debug_x(:, location+1:end));
    else
        sirt = SIRT(fun, model.d, sirtold, 'sample_x', debug_x(:, 1:location), 'debug_x', debug_x(:, location+1:end));
    end
    model.sirt{k} = sirt;
    
    thetasams = eval_irt(sirt, rand(model.d, 1000));
    w = fun(thetasams)./eval_pdf(sirt, thetasams);
    w(isnan(w)) = 0;
    model.fullESS(k) = sum(w)^2/sum(w.^2)/10;
    sam_uni = rand(model.d, model.sam_size);
    sam_old = L(1: model.d, 1: model.d) * eval_irt(sirt, sam_uni) + mu(1: model.d);
    sam_old = That(sam_old);
    sam_old = sam_old(1:model.td+model.dimension, :);

    pdfthold = @(x) nan2zero(pdfth2(x).*...
            eval_pdf_check(sirt, L(1:model.td+ model.dimension, 1:model.td+ model.dimension)\(Thatinv2(x)-mu(1:model.td+ model.dimension))...
            , model.poly) ./model.fu(Thatinv2(x)));
    ndouble = 0;
    % collapse start
    if strcmp(model.coll, 'NA')
        if k == model.steps
            fprintf('ESS at step %d is %.2f /%d \n', k, stepess, 2^ndouble * model.sam_size)
            fprintf('The preESS at step %d is %.2f  \n',k,model.preESS(k))
            fprintf('The fullESS at step %d is %.2f  \n',k,model.fullESS(k))
            fprintf('The collapseESS at step %d is %.2f  \n\n',k, model.collESS(k))
            break
        end
        sam_temp = st_process(sam_old, model);
        sam_new = sam_old;
        sam_new(model.td+1:model.td+model.dimension, :) = sam_temp;
        sam_new(model.td+model.dimension+1:model.d, :) = sam_old(model.td+1:model.td+model.dimension, :);
        w = like(sam_new, model, model.Y(:, k+1));
        ndouble = 0;
        stepess = sum(w)^2/sum(w.^2);
        % if ESS is too small, double the samples to find support of ftt 
        while stepess < model.sam_size/100 && k~=1       
            ndouble = ndouble+1;
            sam_old2 = L(1: model.d, 1: model.d) * eval_irt(sirt, rand(model.d, 2^(ndouble-1) * model.sam_size)) + mu(1: model.d);   
            sam_old2 = That(sam_old2);
            sam_old2 = sam_old2(1:model.td+model.dimension, :);
            sam_temp2 = st_process(sam_old2, model);
            sam_new2 = sam_old2; % later use like function to calculate likelihood. so sam_new is in real space. only x_k part from mu cov can be used.
            sam_new2(model.td + 1 : model.td + model.dimension, :) = sam_temp2;
            sam_new2(model.td + model.dimension+1 : model.d, :) = sam_old2(model.td + 1:model.td + model.dimension, :);
            w2 = like(sam_new2, model, model.Y(:, k+1));
            w = [w,w2]; %#ok<AGROW>
            sam_new = [sam_new, sam_new2]; %#ok<AGROW>
            stepess = sum(w)^2/sum(w.^2);
            if ndouble > 4
                break
            end
        end
        fprintf('precollapse at step %d is %.2f \n', k, stepess)
        model.stepess(k) = stepess;
        sam_new = datasample(sam_new', model.sam_size,'Weights',w);
        sam_new = sam_new';
        
        
        mu_theta = mean(sam_new, 2);
%         L_theta = chol(cov(sam_new') + 1e-8 * eye(model.d))';
        if model.linear == 0
            L_theta = model.amplify*diag(sqrt(var(sam_new, 0, 2)));
        elseif model.linear == 1
            Sigma = (model.c * cov(sam_new')^(-1) + (1- model.c) * eye(size(sam_new,1)))^(-1);
            L_theta = chol(Sigma + 1e-8 * eye(model.d))';
        end
        
        model.L_theta(:, :, k) = L_theta;
        model.mu_theta(:, k) = mu_theta;
 
        Lx = @(x) L_theta*x+mu_theta;
        pdffull = @(x) (nan2zero(pdfth2(rt2*Lx(x)).* ...
            eval_pdf_check(sirt, L(1:model.td+ model.dimension, 1:model.td+ model.dimension)\(Thatinv2(rt2*Lx(x))- mu(1:model.td+ model.dimension)), model.poly)./...
            model.fu(Thatinv2(rt2*Lx(x)))) .*transition(Lx(x), model).*like(Lx(x), model, model.Y(:,k+1)))./model.fu(x);
%         const = pdffull(zeros(model.d, 1))*model.c;
        pdffull = @(x) pdffull(x).^model.c.*model.fu(x);
%         pdffull = @(x) pdffull(x) + const;

        pdfth_app = pdffull;

        debug_x = L_theta\(sam_new-mu_theta);
        location = floor(size(debug_x,2)/2);


        fprintf('\n collapse for Tth at time %d\n', k)
        if k == 2
            reirt = SIRT(pdffull, model.d, model.polycoll, model.optcoll, 'sample_x', debug_x(:, 1:location), 'debug_x', debug_x(:, location+1:end));
        else
            reirt = SIRT(pdffull, model.d, reirt,  'sample_x', debug_x(:, 1:location), 'debug_x', debug_x(:, location+1:end));
        end
        model.reirt{k} = reirt;

        That = @(u) L_theta * eval_irt(reirt, model.R(u)) + mu_theta;
        Thatinv = @(x) model.Rinv(eval_rt(reirt, L_theta\(x-mu_theta)));
        Thatinv2 = @(x) model.Rinv(eval_rt(reirt, L_theta(1:model.td+model.dimension, 1:model.td+model.dimension)\(x-mu_theta(1:model.td+model.dimension))));
        pdfth = @(x) nan2zero(eval_pdf_check(reirt, L_theta\(x-mu_theta), model.polycoll));
        pdfth2 = @(x) nan2zero(eval_pdf_check(reirt, L_theta(1:model.td+model.dimension, 1:model.td+model.dimension)\(x-mu_theta(1:model.td+model.dimension)), model.polycoll));
        
        thetasams = eval_irt(reirt, rand(model.d, 1000));
        w = pdfth_app(thetasams)./eval_pdf(reirt, thetasams);
        model.collESS(k) = sum(w)^2/sum(w.^2)/10;
       
        sam_L = Thatinv(sam_new(1:model.d,:));
    else
        error('no collapse precondition option')
    end
    % collapse end
    
    
    fprintf('ESS at step %d is %.2f /%d \n', k, stepess, 2^ndouble * model.sam_size)
    fprintf('The preESS at step %d is %.2f  \n',k,model.preESS(k))
    fprintf('The fullESS at step %d is %.2f  \n',k,model.fullESS(k))
    fprintf('The collapseESS at step %d is %.2f  \n\n',k, model.collESS(k))
    model.time_filter(k) = toc;
    fprintf('this iteration took %.2f s \n', model.time_filter(k))
end                         
end