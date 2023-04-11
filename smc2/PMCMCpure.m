function [model, PF, logweight_t, para, ar] = PMCMCpure(model, PF, logweight_t, T)
%{
the PMCMC step in SMC^2
T is the current time of states
for each parameter, run an independent MC. Pick the first value after
burnin
Apply the adapted version of MCMC
%}
NN = 5000;
model.Nx = 2*model.Nx;
PF_star = cell(model.sam_size, 1); % each entry contains filter states in an PF for the para, of dimension d*N


para = zeros(model.td*model.sam_size, model.burnin + NN);
reverseStr = [];
ar = [];
for num = 1: model.burnin + NN % each loop is a MCMC
    % propose new samples
    sigma = chol(cov(model.para') + 1e-5*eye(model.td));
    para_star = model.para + sigma * randn(model.td, model.sam_size);
   
    logweight_star = zeros(model.sam_size, 1);
    for k = 1: model.sam_size
       % associate PF to each proposed parameters
       reverseStr = displayprogress(100*((num-1)*model.sam_size + k)/(model.burnin + NN)/model.sam_size, reverseStr);
       theta = para_star(:, k);
       X_old = samx0(theta, model.Nx, model); % dimension d*N      
       for time = 1:T
           X_new = st_process([repmat(theta, 1, model.Nx);X_old], model);
           weight_new = like([repmat(theta, 1, model.Nx);X_new], model, model.Y(:, time));
           logweight_star(k, 1) = logweight_star(k, 1) + log(mean(weight_new));
           X_new = datasample(X_new, model.Nx, 2, 'Weight', weight_new);
           X_old = X_new;
       end  
       PF_star{k} = X_new;
       % acceptance and rejection
       logacceptprob = logweight_star(k, 1) + log(parapriorpdf(theta, model)) ...
           - logweight_t(k, 1) - log(parapriorpdf(model.para(:, k), model));
%         logacceptprob
       
       if log(rand(1)) < logacceptprob
          model.para(:, k) = para_star(:, k);
          PF{k} = PF_star{k};
          logweight_t(k, 1) = logweight_star(k, 1);
          ar = [ar,1];
       else
          ar = [ar,0];
       end
    end
   para(:, num) = reshape(model.para, [], 1);
end
logweight_t = zeros(model.sam_size, 1);
end

function reverseStr = displayprogress(perc,reverseStr)
msg = sprintf('%3.1f', perc);
fprintf([reverseStr, msg, '%%']);
reverseStr = repmat(sprintf('\b'), 1, length(msg)+1);
end
