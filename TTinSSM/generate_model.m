function model = generate_model(model)
%  GENERATE_MODEL Given the model parameters, generate data of the model, together with necessary 
% functions needed.
% 
% model.name : kalman, lorenz, SV, PP
% 
% model.dimension : states dimension
% 
% model.td : theta dimension
% 
% model.X : states
% 
% model.Y : observations
% 
% model.steps : steps for data generation
% 
% model.initg : prior
% 
% model.stepg : update
% 
% model.pdf : theoretical pdf
% 
% model.cmat : transformation matrix from [theta, t+1, t] to [theta, t]
% 
% model.ftt : store the ftt in every step
% 
% model.mu_all , model.L_all : record all the transformations needed
% 
% model.d : build ftt in d dimension
% 
% 
% input all funcitons needed to data structure model.
% generate states and observations
% initialise all the tensors
% initialise all the tensors
model.Y = zeros(size(model.C,1), model.steps); 
model.X = zeros(model.dimension, model.steps+1); 

model.sirt = cell(model.steps, 1);
model.reirt = cell(model.steps, 1); % store the re-approximated pdf of theta
model.d = model.td + 2*model.dimension;
if strcmp(model.pre, 'DIRT') && strcmp(model.coll, 'adapted')
    model.colltime = 1;
    model.collind = 1:model.steps; 
    model.collind = max(model.collind - model.colltime, 1); % store the collapse precondition index for each iteration
    model.collirt = cell(model.steps, 1); % store the direct irt form u_\theta to \theta, used for preconditioning in collapse step
end
model.cmat = [eye(model.td, model.td + 2*model.dimension) ; [zeros(model.dimension, model.td+model.dimension), eye(model.dimension)]];
model.xmat = [zeros(model.dimension, model.td),eye(model.dimension)];
model.mu_all = zeros(model.td + 2*model.dimension, model.steps);
model.L_all = zeros(model.td + 2*model.dimension, model.td + 2*model.dimension, model.steps);

if strcmp(model.pre, 'DIRT') 
    model.mu_theta = zeros(model.d, model.steps);
    model.L_theta = zeros(model.d, model.d, model.steps);
    tail = model.tail;
    model.fu = @(x) tg_pdf(x, tail);
    model.R = @(x) tg_cdf(x, tail);
    model.Rinv = @(x) tg_inv(x, tail);
elseif strcmp(model.pre, 'SSM') 
    model.mu_theta = zeros(model.td, model.steps);
    model.L_theta = zeros(model.td, model.td, model.steps);
    tail = model.tail;
    model.fu = @(x) tg_pdf(x, tail);
    model.R = @(x) tg_cdf(x, tail);
    model.Rinv = @(x) tg_inv(x, tail);
end


% time record
model.ess = zeros(model.steps, 1);
model.stepess = zeros(model.steps, 1);
model.preESS = zeros(model.steps, 1);
model.fullESS = zeros(model.steps, 1);
model.collESS = zeros(model.steps, 1);
model.time_filter = zeros(model.steps, 1);
model.time_sample = zeros(model.steps, 1);


if  strcmp(model.pre, 'DIRT') && strcmp(model.coll, 'adapted')
    model.mu_coll = zeros(model.td, model.steps);
    model.L_coll = zeros(model.td, model.td, model.steps);
end

% generate states and observations
sams = priorsam(model);
model.X(:,1) = sams(model.td+1 : model.td+model.dimension, 1);

for k = 1: model.steps
    model.X(:, k+1) = st_process([model.theta_transformed; model.X(:,k)], model); 
    model.Y(:, k) = ob_process(model.X(:, k+1), model);
end
end