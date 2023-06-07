%% smoothing procedure based on the TT filtering result

% clearvars -except rootref
% clf
file = 'test';
root = pwd;
load([root '/data/' file '.mat'])

timeindex = 1:model.steps;  % the time index to compute ESS
NN = 40; % number of batches used to generate the paths
model.N = 5e3; % sample size to compute ESS 5e3
addpath([rootref '/' model.name])


%% ESS

ESS = zeros(NN, model.steps);
time_sample = zeros(NN, model.steps);


for t = timeindex
   if strcmp(model.pre, 'DIRT')
       stats = generate_samples_correctionfull(model,t, NN);
   elseif strcmp(model.pre, 'NL')
       stats = generate_samples_nlog(model,t, NN);
%    elseif strcmp(model.pre, 'hybrid')
%        stats = generate_samples_hybrid(model,t, NN);
   elseif strcmp(model.pre, 'SSM')
       stats = generate_samples_dirtlog(model,t, NN);
   end
   fprintf(' %d over %d has been finished\n', t, timeindex(end))
   ESS(:, t) = stats.ess;
   time_sample(1, t) = stats.time_sample;
end

model.ess = mean(ESS(:, timeindex));
figure(2)
errorbar(timeindex, mean(ESS(:, timeindex)), std(ESS(:, timeindex)));
title('ESS')






save([root '/data/' file])


%%
rmpath([rootref '/' model.name])