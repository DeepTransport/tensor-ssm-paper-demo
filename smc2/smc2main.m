%% Sequential transport for kalman test
%{
this is the body function for SMC2
All the parameters are defined in the first section
the output are parameter smaples
%}
clearvars -except rootref
clc
%% Model parameters
model.name = 'kalman';
tic
model.steps = 50;   % the steps for sequential inference
                   % for cts model, the total time is model.dt * steps 
model.sam_size = 1e4; % the sample size for parameters
model.Nx = 1e2; % the initial sample size for X-PF
PF = cell(model.sam_size, 1); % each entry contains filter states in an PF for the para, of dimension d*N
logweight_t = zeros(model.sam_size, 1); % the weights associated with each parameter

model.gamma = 0.5; % the criteria to trigger PMCMC
model.burnin = 1;

file = 'kalman'; % the file folder storing results
model.file = file;


%% complete model
mkdir('results', file)
addpath([rootref '\' model.name])

rng(1)
model = truetheta(model);   % input parameters given in truetheta to the model
model = generate_model(model);   % generate states and observations using thetas given. 
rng('shuffle')


% load data50kalman % data contains 50 step kalmanfilter with theta [0.8,0.5]
% model.C = C;
% model.X = X(:, 1:model.steps+1);
% model.Y = Y(:, 1:model.steps);



% load data1000sv 
% model.C = C;
% model.X = X;
% model.Y = Y;

%% initialization
sams = priorsam(model);
model.paraall = zeros(model.td, model.sam_size, model.steps+1);
model.para = sams(1:model.td, 1:model.sam_size); % this is the parameter posterior samples. each col is a sample
model.paraall(:, :, 1) = model.para;
model.weight = ones(model.sam_size, model.steps+1);

% time 0 initialization
for k = 1: model.sam_size
   X0 = samx0(model.para(:, k), model.Nx, model); % dimension d*N
   PF{k} = X0;
end

for t = 1:model.steps % the smc for parameters
    % forward PF for each parameter sample
    fprintf('time is %d now\n', t);
    for k = 1: model.sam_size
       theta = model.para(:, k);
       X_old = PF{k};
       
       X_new = st_process([repmat(theta, 1, model.Nx); repmat(X_old, 1, model.Nx/size(X_old, 2))], model);
       weight_new = like([repmat(theta, 1, model.Nx);X_new], model, model.Y(:, t));
       logweight_t(k, 1) = logweight_t(k, 1) + log(mean(weight_new));
       X_new = datasample(X_new, model.Nx, 2, 'Weight', weight_new);
       
       PF{k} = X_new;
    end
    
    % trigger PMCMC when ESS falls below criteria
    weight = exp(logweight_t);
    if sum(weight)^2/sum(weight.^2)/model.sam_size < model.gamma
        fprintf('the MCMC step begins at %d, the number of particles is doubled from %d\n', t, model.Nx)
        [model, PF, logweight_t] = PMCMC(model, PF, logweight_t, t);    
    end
    model.paraall(:,:,t+1) = model.para;
    model.weight(:, t+1) = exp(logweight_t);
end
model.time = toc;
smc2 = model;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PMCMC
model.sam_size = 60;
model.burnin = 300;
sams = priorsam(model);
model.paraall = zeros(model.td, model.sam_size, model.steps+1);
model.para = sams(1:model.td, 1:model.sam_size); % this is the parameter posterior samples. each col is a sample
model.paraall(:, :, 1) = model.para;
model.weight = ones(model.sam_size, model.steps+1);

% initialization
for k = 1: model.sam_size
   X0 = samx0(model.para(:, k), model.Nx, model); % dimension d*N
   PF{k} = X0;
end

for t = 1:model.steps % the smc for parameters
    % forward PF for each parameter sample
    fprintf('time is %d now\n', t);
    for k = 1: model.sam_size
       theta = model.para(:, k);
       X_old = PF{k};
       
       X_new = st_process([repmat(theta, 1, model.Nx); repmat(X_old, 1, model.Nx/size(X_old, 2))], model);
       weight_new = like([repmat(theta, 1, model.Nx);X_new], model, model.Y(:, t));
       logweight_t(k, 1) = logweight_t(k, 1) + log(mean(weight_new));
       X_new = datasample(X_new, model.Nx, 2, 'Weight', weight_new);
       
       PF{k} = X_new;
    end
end


[model, PF, logweight_t, paramcmc, ar] = PMCMCpure(model, PF, logweight_t, model.steps);
model.time = toc;
%% plot for kalman filter
N = 100;
index = linspace(0.4, 1 ,N);
index2 = index + 1e-2*eye(1,N) - 1e-2* flip(eye(1,N));
[X, Y] = meshgrid(index,index);
[X_index2, Y_index2] = meshgrid(index2,index2);
vgrid = [X(:),Y(:)]';


para = normcdf(model.paraall(:,:, end))*0.6+0.4;
[f, xi] = ksdensity(para', vgrid', 'Weights', model.weight(:, end));

figure
% contour(reshape(f,N,N))
mesh(X,Y,reshape(f,N,N))
title('SMC2')


p_thm = reshape(theta41(vgrid, model.C, model.Y(:, 1:model.steps)), N, N);
p_thm = p_thm./(nansum(nansum(p_thm))* (index(2)-index(1))^2);
figure
% contour(p_thm)
mesh(X,Y,p_thm)
title('theory')
    


%%
root = pwd;
save([root '\results\' file '\' file],'model');
rmpath([rootref '\' model.name])

