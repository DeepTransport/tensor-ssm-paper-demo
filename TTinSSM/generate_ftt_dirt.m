 function model = generate_ftt_dirt(model)
%  GENERATE_FTT Given the model, generate ftt in every step.

% define some transition matrix
rt = sparse([eye(model.td, model.td), zeros(model.td, model.dimension * 2)]); % transition matrix to get theta
r1 = sparse([zeros(model.dimension, model.td), eye(model.dimension),zeros(model.dimension,model.dimension)]); % transition matrix to get x1
rt2 = sparse(model.td+model.dimension, model.d);
rt2(1:model.td, 1:model.td) = eye(model.td); % transition matrix to get theta and x2
rt2(model.td+1:model.td+model.dimension, model.d-model.dimension+1:model.d) = eye(model.dimension);  % transition matrix to get theta and x2


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

L = chol(cov(sam_new') + 1e-8 * eye(size(sam_new,1)))'; % linear transformation
mu = mean(sam_new, 2);

% store block matrices
model.L_all(:, :, 1) = L;
model.mu_all(:, 1) = mu;
% regularization ends

fun = @(x) fun(model.cmat*(L*x+mu)).*transition((L*x+mu), model).*...
        like((L*x+mu), model, model.Y(:,1)); % update
sirt = SIRT(fun, model.d, model.poly, model.optinit);
model.sirt{1} = sirt;

% build T_{1}^x
Tx = @(u) model.xmat * (L(1: model.td+model.dimension, 1: model.td+model.dimension) * eval_irt(sirt, model.R(u)) + mu(1: model.td+model.dimension));

% new samples
sam_uni = rand(model.d - model.dimension, model.sam_size);
sam_u = model.Rinv(sam_uni); % this will be the input samples to next iteration(used to calculate L) 
sam_old = L(1: model.td+model.dimension, 1: model.td+model.dimension) * eval_irt(sirt, sam_uni) + mu(1: model.td+model.dimension);

% build T_{1}^\theta
Tth = @(u) L(1: model.td, 1: model.td) * eval_irt(sirt, model.R(u)) + mu(1: model.td);
Tthinv = @(x) model.Rinv(eval_rt(sirt, L(1: model.td, 1: model.td)\(x-mu(1: model.td))));
pdfth = @(x) eval_pdf_check(sirt, L(1: model.td, 1: model.td) \ (x - mu(1:model.td)), model.poly);
% the following three functions are the base transport of the future Tth
% functions
Tthbs = Tth;
Tthinvbs = Tthinv;
pdfthbs = pdfth;
model.time_filter(1) = toc;

if strcmp(model.coll, 'adapted')
    model.collirt{1} = sirt;
    model.mu_coll(:, 1) = mu(1:model.td);
    model.L_coll(:,:, 1) = L(1:model.td, 1:model.td);
end

for k = 2 : model.steps
    tic
    % regularization begins
    sam_temp = st_process(sam_old, model);
    sam_new = sam_old; % later use like function to calculate likelihood. so sam_new is in real space. only x_k part from mu cov can be used.
    sam_new(model.td + 1 : model.td + model.dimension, :) = sam_temp;
    sam_new(model.td + model.dimension+1 : model.d, :) = sam_old(model.td + 1:model.td + model.dimension, :);
    w = like(sam_new, model, model.Y(:, k));
    ndouble = 0;
    stepess = sum(w)^2/sum(w.^2);
    % if ESS is too small, double the samples to find support of ftt 
    while stepess < model.sam_size/100 && k~=2        
        ndouble = ndouble+1;
        sam_uni2 = rand(model.d - model.dimension, 2^(ndouble-1) * model.sam_size);
        sam_u2 = model.Rinv(sam_uni2); % this will be the input samples to next iteration(used to calculate L) 
        sam_old2 = L(1: model.td+model.dimension, 1: model.td+model.dimension) * eval_irt(sirt, sam_uni2) + mu(1: model.td+model.dimension);
        sam_old2(1: model.td, :) = Tth_pre(sam_old2(1: model.td, :));        
        sam_temp2 = st_process(sam_old2, model);
        sam_new2 = sam_old2; % later use like function to calculate likelihood. so sam_new is in real space. only x_k part from mu cov can be used.
        sam_new2(model.td + 1 : model.td + model.dimension, :) = sam_temp2;
        sam_new2(model.td + model.dimension+1 : model.d, :) = sam_old2(model.td + 1:model.td + model.dimension, :);
        w2 = like(sam_new2, model, model.Y(:, k));
        w = [w,w2]; %#ok<AGROW>
        sam_temp = [sam_temp, sam_temp2]; %#ok<AGROW>
        sam_u = [sam_u, sam_u2]; %#ok<AGROW>
        stepess = sum(w)^2/sum(w.^2);
        if ndouble > 4
            break
        end
    end
    
    model.stepess(k) = stepess/2^(ndouble-1)/model.sam_size*100; % normalized to 100
    
    sam_u_re = [sam_u; sam_u(model.td+1:model.td+model.dimension, :)];
    sam_u_re(model.td+1:model.td+model.dimension, :) = sam_temp;
    
    sam_new = datasample(sam_u_re', model.sam_size,'Weights',w);
    sam_new = sam_new'; % the new sam_new is d * samsize 

    
    sam_new(:, any(isnan(sam_new))) = [];
    L = chol(cov(sam_new') + 1e-8)';
    mu = mean(sam_new, 2);
    model.L_all(:, :, k) = L;
    model.mu_all(:, k) = mu;
    debug_x = L\(sam_new - mu);
    location = floor(size(debug_x,2)/2);
    % regularization ends
    
    % computation of preESS
    Z = randn(model.d, 1000);
    LZ = L*Z + mu;
    X = LZ;
    X(1:model.td, :) = Tth(X(1:model.td, :));
    X( model.td+model.dimension+1 : model.d, :) = Tx(LZ( [1:model.td, model.td+model.dimension+1:model.d] ,:));
    w = transition(X, model).*like(X, model, model.Y(:, k)).* ...
        model.fu(LZ([1:model.td, model.td+model.dimension+1 : model.d],:))./model.fu(Z);
    model.preESS(k) = sum(w)^2/sum(w.^2)/10;

    fprintf('this is ftt build up at time %d\n', k)
    if k == 2
        fun = @(u) transition([Tth(rt*(L*u+mu));r1*(L*u+mu);Tx(rt2*(L*u+mu))], model).*...
         like([Tth(rt*(L*u+mu));r1*(L*u+mu);Tx(rt2*(L*u+mu))], model, model.Y(:,k)).*model.fu(rt2*(L*u+mu));   
    else
        fun = @(u) transition([Tth(rt*(L*u+mu));r1*(L*u+mu);Tx(rt2*(L*u+mu))], model).*...
         like([Tth(rt*(L*u+mu));r1*(L*u+mu);Tx(rt2*(L*u+mu))], model, model.Y(:,k)).*model.fu(rt2*(L*u+mu)) .* ...
         pthold(Tth(rt*(L*u+mu)))./zero2inf(pdfth(Tth(rt*(L*u+mu)))); % update, third line is correction
    end
    if k == 2 
        sirt = SIRT(fun, model.d, model.poly, model.opt, 'sample_x', debug_x(:, 1:location), 'debug_x', debug_x(:, location+1:end));
    else
        sirt = SIRT(fun, model.d, sirt, 'sample_x', debug_x(:, 1:location), 'debug_x',debug_x(:, location+1:end));
    end
    model.sirt{k} = sirt;
    
    % computation of fullESS
    Z = randn(model.d, 1000);
    RZ = model.R(Z);
    [SRZ, SRZf] = eval_irt(sirt, RZ);
    LSRZ = L*SRZ + mu;
    X = LSRZ;
    X(1:model.td, :) = Tth(LSRZ(1:model.td, :));
    X(model.td+model.dimension+1 : model.d, :) = Tx(LSRZ( [1:model.td, model.td+model.dimension+1:model.d] ,:));
    w = transition(X, model).*like(X, model, model.Y(:, k)).* ...
        model.fu(LSRZ([1:model.td, model.td+model.dimension+1 : model.d], :))./SRZf; %eval_pdf(sirt, SRZ);
    w(isnan(w)) = 0;
    model.fullESS(k) = sum(w)^2/sum(w.^2)/10;

%     sam_theta_last = sam_old(1:model.td, :); % save the theta values from last iteration, in case the collapse pre is needed
    % generate samples and map to space X
    sam_uni = rand(model.d - model.dimension, model.sam_size);
    sam_u = model.Rinv(sam_uni); % this will be the input samples to next iteration(used to calculate L) 
    sam_old = L(1: model.td+model.dimension, 1: model.td+model.dimension) * eval_irt(sirt, sam_uni) + mu(1: model.td+model.dimension);
    sam_old(1: model.td, :) = Tth(sam_old(1: model.td, :));
    Tx = @(u) model.xmat * (L(1: model.td+model.dimension, 1: model.td+model.dimension) * eval_irt(sirt, model.R(u)) + mu(1: model.td+model.dimension));
    Tthold = @(u) Tth(L(1: model.td, 1: model.td) * eval_irt(sirt, model.R(u)) + mu(1: model.td));
    if k == 2
        pthold = @(x) nan2zero(pdfth(x).*...
                eval_pdf_check(sirt, L(1:model.td, 1:model.td)\(Tthinv(x)-mu(1:model.td)), model.poly)...
                     ./model.fu(Tthinv(x)));
    else
        pthold = @(x) nan2zero(pdfth(x).*...
                eval_pdf_check(sirt, L(1:model.td, 1:model.td)\(Tthinv(x)-mu(1:model.td)), model.poly)...
                     ./model.fu(Tthinv(x)));
    end
    
    % collapse start
    if strcmp(model.coll, 'NA')
        mu_theta = mean(sam_old(1:model.td, :), 2);
        L_theta = chol(cov(sam_old(1:model.td, :)') + 1e-8 * eye(model.td))';
        model.L_theta(:, :, k) = L_theta;
        model.mu_theta(:, k) = mu_theta;
        if k == 2
            pdfth = @(x) nan2zero(pdfth(L_theta*x+mu_theta).*...
                eval_pdf_check(sirt, L(1:model.td, 1:model.td)\(Tthinv(L_theta*x+mu_theta)-mu(1:model.td)), model.poly)...
                     ./model.fu(Tthinv(L_theta*x+mu_theta)));
        else
            pdfth = @(x) nan2zero(eval_pdf(model.reirt{k-1}, model.L_theta(:,:,k-1)\(L_theta*x+mu_theta - model.mu_theta(:, k-1))).*...
                eval_pdf_check(sirt, L(1:model.td, 1:model.td)\(Tthinv(L_theta*x+mu_theta)-mu(1:model.td)), model.poly)...
                     ./model.fu(Tthinv(L_theta*x+mu_theta)));
        end
        pdfthold = pdfth;

        debug_x = L_theta\(sam_old(1: model.td, :)-mu_theta);

        fprintf('\n collapse for Tth at time %d\n', k)
        if k == 2
            reirt = SIRT(pdfth, model.td, model.polycoll, model.optcoll, 'sample_x', debug_x(:, 1:model.sam_size/2), 'debug_x', debug_x(:, model.sam_size/2:model.sam_size));
        else
            reirt = SIRT(pdfth, model.td, reirt,  'sample_x', debug_x(:, 1:model.sam_size/2), 'debug_x', debug_x(:, model.sam_size/2:model.sam_size));
        end
        model.reirt{k} = reirt;

        Tth_pre = Tth;
        Tth = @(u) L_theta * eval_irt(reirt, model.R(u)) + mu_theta;
        Tx = @(u) model.xmat * (L(1: model.td+model.dimension, 1: model.td+model.dimension) * eval_irt(sirt, model.R(u)) + mu(1: model.td+model.dimension));
        Tthinv = @(x) model.Rinv(eval_rt(reirt, L_theta\(x-mu_theta)));
        pdfth = @(x) nan2zero(eval_pdf_check(reirt, L_theta\(x-mu_theta), model.polycoll));
        
        thetasams = eval_irt(reirt, rand(model.td, 1000));
        w = pdfthold(thetasams)./eval_pdf(reirt, thetasams);
        model.collESS(k) = sum(w)^2/sum(w.^2)/10;
        
    elseif strcmp(model.coll, 'time1')
        uth = Tthinvbs(sam_old(1:model.td, :));
        uth(:, any(isnan(uth))) = [];
        
        mu_theta = mean(uth, 2);
        L_theta = chol(cov(uth') + 1e-8 * eye(model.td))';
        model.L_theta(:, :, k) = L_theta;
        model.mu_theta(:, k) = mu_theta;

        pdfth_u = @(u) nan2zero(pdfth(Tthbs(L_theta*u+mu_theta)).*...
                eval_pdf_check(sirt, L(1:model.td, 1:model.td)\(Tthinv(Tthbs(L_theta*u+mu_theta))-mu(1:model.td)), model.poly)...
                     ./model.fu(Tthinv(Tthbs(L_theta*u+mu_theta))).* model.fu(L_theta*u+mu_theta)./ pdfthbs(Tthbs(L_theta*u+mu_theta)));

        debug_x = L_theta\(uth-mu_theta);
        location = floor(size(debug_x,2)/2);

        fprintf('\n collapse for Tth at time %d\n', k)
        if k == 2
            reirt = SIRT(pdfth_u, model.td, model.polycoll, model.optcoll, 'sample_x', debug_x(:, 1:location), 'debug_x', debug_x(:, location+1:end));
        else
            reirt = SIRT(pdfth_u, model.td, reirt, 'sample_x', debug_x(:, 1:location), 'debug_x', debug_x(:, location+1:end));
        end
        model.reirt{k} = reirt;

        Tth_pre = Tth;
        Tth = @(u) Tthbs(L_theta * eval_irt(reirt, model.R(u)) + mu_theta);
        Tthinv = @(x) model.Rinv(eval_rt_check(reirt, L_theta\(Tthinvbs(x)-mu_theta), model.poly));
        pdfth = @(x) nan2zero(eval_pdf_check(reirt, L_theta\(Tthinvbs(x)-mu_theta), model.polycoll).* pdfthbs(x) ./...
        model.fu(Tthinvbs(x)));
    
        thetasams = eval_irt(reirt, rand(model.td, 1000));
        w = pdfth_u(thetasams)./eval_pdf(reirt, thetasams);
        model.collESS(k) = nansum(w)^2/nansum(w.^2)/10;
        
    elseif strcmp(model.coll, 'adapted')
        % choose base function
        collind = model.collind(k);
        Tthbs = @(u) model.L_coll(1: model.td, 1: model.td, collind) * eval_irt(model.collirt{collind}, model.R(u)) + model.mu_coll(1: model.td, collind);
        Tthinvbs = @(x) model.Rinv(eval_rt(model.collirt{collind}, model.L_coll(1: model.td, 1: model.td, collind)\(x-model.mu_coll(1: model.td, collind))));
        pdfthbs = @(x) eval_pdf_check(model.collirt{collind}, model.L_coll(1: model.td, 1: model.td, collind) \ (x - model.mu_coll(1: model.td, collind)), model.poly);
        uth = Tthinvbs(sam_old(1:model.td, :));
        uth(:, any(isnan(uth))) = [];
        mu_theta = mean(uth, 2);
        L_theta = chol(cov(uth') + 1e-8 * eye(model.td))';
        model.L_theta(:, :, k) = L_theta;
        model.mu_theta(:, k) = mu_theta;

        pdfth_u = @(u) nan2zero(pdfth(Tthbs(L_theta*u+mu_theta)).*...
                eval_pdf_check(sirt, L(1:model.td, 1:model.td)\(Tthinv(Tthbs(L_theta*u+mu_theta))-mu(1:model.td)), model.poly)...
                     ./model.fu(Tthinv(Tthbs(L_theta*u+mu_theta))).* model.fu(L_theta*u+mu_theta)./ pdfthbs(Tthbs(L_theta*u+mu_theta)));

        debug_x = L_theta\(uth-mu_theta);

        fprintf('\n collapse for Tth at time %d\n', k)
        if k == 2
            reirt = SIRT(pdfth_u, model.td, model.polycoll, model.optcoll, 'sample_x', debug_x(:, 1:floor(size(debug_x, 2)/2)), 'debug_x', debug_x(:, floor(size(debug_x, 2)/2):size(debug_x, 2)));
        else
            try	
            reirt = SIRT(pdfth_u, model.td, reirt, 'sample_x', debug_x(:, 1:floor(size(debug_x, 2)/2)), 'debug_x', debug_x(:, floor(size(debug_x, 2)/2):size(debug_x, 2)-1));
            catch
            save('errorfile', 'pdfth_u', 'reirt', 'debug_x')
            error('pdfth_u')
            end
        end
        model.reirt{k} = reirt;
        
        thetasams = eval_irt(reirt, rand(model.td, 1000));
        w = pdfth_u(thetasams)./eval_pdf(reirt, thetasams);
        model.collESS(k) = sum(w)^2/sum(w.^2)/10;
        
        Tth_pre = Tth;
        Tth = @(u) Tthbs(L_theta * eval_irt(reirt, model.R(u)) + mu_theta);
        Tthinv = @(x) model.Rinv(eval_rt_check(reirt, L_theta\(Tthinvbs(x)-mu_theta), model.poly));
        pdfth = @(x) nan2zero(eval_pdf_check(reirt, L_theta\(Tthinvbs(x)-mu_theta), model.polycoll).* pdfthbs(x) ./...
        model.fu(Tthinvbs(x)));
        % here the renewal of collapse preconditioner begins

        
        L_coll = chol(cov(sam_old(1:model.td, :)') + 1e-8 * eye(model.td))';
        mu_coll = mean(sam_old(1:model.td, :), 2);
        debug_x = L_coll\(sam_old(1:model.td, :)-mu_coll);
        pdf_coll = @(x) nan2zero(pdfth(L_coll*x + mu_coll));
%         pdf_coll = @(x) nan2zero(pdfth(L_coll*x + mu_coll).*...
%                 eval_pdf_check(sirt, L(1:model.td, 1:model.td)\(Tthinv(L_coll*x + mu_coll)-mu(1:model.td)), model.poly)...
%                      ./model.fu(Tthinv(L_coll*x + mu_coll)));
        if k == 2
            collirt = SIRT(pdf_coll, model.td, model.polycoll, model.opt, 'sample_x', debug_x(:, 1:model.sam_size/2), 'debug_x', debug_x(:, model.sam_size/2:model.sam_size));
        else
            collirt = SIRT(pdf_coll, model.td, model.collirt{k-1}, 'sample_x', debug_x(:, 1:model.sam_size/2), 'debug_x', debug_x(:, model.sam_size/2:model.sam_size));
        end

        model.collirt{k} = collirt;
        model.mu_coll(:, k) = mu_coll;
        model.L_coll(:, :, k) = L_coll;
    else
        error('no collapse precondition option')
    end
    % collapse end
    
    try
        save([root '\ppresult\' file '\modelstep',num2str(k)], 'model')
        save([root '\ppresult\' file '\tempstep',num2str(k)])
    catch
    end
    
    fprintf('ESS at step %d is %.2f /%d \n', k, stepess, 2^ndouble * model.sam_size)
    fprintf('The preESS at step %d is %.2f  \n',k,model.preESS(k))
    fprintf('The fullESS at step %d is %.2f  \n',k,model.fullESS(k))
    fprintf('The collapseESS at step %d is %.2f  \n\n',k, model.collESS(k))
    model.time_filter(k) = toc;
    fprintf('this iteration took %.2f s \n', model.time_filter(k))
end                         
end