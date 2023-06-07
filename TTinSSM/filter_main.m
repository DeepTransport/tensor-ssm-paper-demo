%% Sequential inference in SSM
%{
this is the body function for tensor-train approach in state space models
Including the basic version which truncates the target function in a bounded domain
Two precondtion versions are provided: NL for nonlinear and DIRT for IRT precondition
In DIRT method, model.c controls the power of the tempering measure
%}

% Model parameters
model.name = 'kalman';   % model names: 'kalman','pp', 'sv'
model.pre = 'NL'; % three different options: NL and DIRT
model.steps = 20;   % the steps for sequential inference
                   % for cts model, the total time is model.dt * steps 
model.sam_size = 1e4; % the sample size to estimate precondition in filter

file = 'test';
model.file = file;



if strcmp(model.pre, 'DIRT') || strcmp(model.pre, 'SSM')
    model.tail = 4;
    model.coll = 'NA'; % in DIRT version, collapse is needed. NA doesn't apply precondition
                     % time1 use parameter posterior at time 1 as preconditioner
                     % adapted use adapted parameter posterior as preconditioner
    model.tail = 4;
    model.c = 0.4; % preconditioning power
    model.amplify = 1.25;
    model.linear = 1;
    model.predim = 'app';
elseif strcmp(model.pre, 'NL')
    model.trans = 'powertail';   % options: powertail and truncated
    model.quantile = 0.002;   % quantiles to set L in powertail method
    model.sdn = 3;   % number of sd tp set L in truncated method
    model.pow = 1;   % choose which powertail function
end


%% complete model
mkdir('data')
addpath([rootref '/' model.name])

rng(1)
model = truetheta(model);   % input parameters given in truetheta to the model
model = generate_model(model);   % generate states and observations using thetas given. states are observations are generated up to time 100
rng('shuffle')

% load data50kalman % data contains 50 step kalmanfilter with theta [0.8,0.5]
% model.C = C;
% model.X = X(:, 1:model.steps+1);
% model.Y = Y(:, 1:model.steps);

%% ftt
poly1 = Legendre(40, [-4, 4]);
poly2 = Lagrange1(50, [-4, 4]);
poly3 = Lagrangep(4, 8, [0, 1], 'ghost_size', 1E-2, 'bc', 'Neumann');
poly4 = Fourier(40, [-4,4]);
poly5 = Lagrange1(60, [-4, 4], 'ghost_size', 1E-3, 'bc', 'Neumann');
poly6 = Lagrangep(2, 20, [-4, 4], 'ghost_size', 1E-3, 'bc', 'Neumann');
poly7 = Lagrange1(80, [-4, 4], 'ghost_size', 1E-3, 'bc', 'Neumann');
poly8 = Lagrangep(4, 8, [-4, 4], 'ghost_size', 1E-3, 'bc', 'Neumann');



opt1 = FTToption('tt_method', 'random', 'max_als', 6, 'als_tol', 1E-5, 'local_tol', 1E-5, 'kick_rank', 5, 'init_rank', 5, 'max_rank', 15, 'sqrt_flag',true);
opt2 = FTToption('tt_method', 'amen', 'max_als', 6, 'als_tol', 1E-5, 'local_tol', 1E-5, 'kick_rank', 5, 'init_rank', 5, 'max_rank', 15, 'sqrt_flag',true);
opt3 = FTToption('tt_method', 'amen', 'max_als', 6, 'als_tol', 1E-5, 'local_tol', 1E-5, 'kick_rank', 5, 'init_rank', 5, 'max_rank', 15, 'sqrt_flag',true);
opt4 = FTToption('tt_method', 'amen', 'max_als', 6, 'als_tol', 1E-5, 'local_tol', 1E-5, 'kick_rank', 5, 'init_rank', 5, 'max_rank', 15, 'sqrt_flag',true);


model.poly = poly3;
model.polycoll = poly6;
model.opt = opt1;
model.optinit = opt2;
model.optcoll = opt3;
model.optpre = opt4;


if strcmp(model.pre, 'DIRT') || strcmp(model.pre, 'SSM')
    model.poly = poly8;
    model.polycoll = poly8;
    model.opt = opt1;
    model.optinit = opt2;
    model.optcoll = opt3;
    model.optpre = opt4;
elseif strcmp(model.pre, 'NL')
    model.poly = poly3;
    model.opt = opt1;
    model.optinit = opt2;
end

%%
% check the poly and option
if strcmp(model.pre, 'DIRT') || strcmp(model.pre, 'SSM')
    if model.poly.domain(1) ~= -model.tail || model.poly.domain(2) ~= model.tail
        error('model polynomial domain is wrong')
    elseif model.polycoll.domain(1) ~= -model.tail || model.polycoll.domain(2) ~= model.tail
        error('model collapse polynomial domain is wrong')
    end
end

if strcmp(model.pre, 'DIRT')
    if strcmp(model.predim, 'app')
    model = generate_ftt_correctionapp(model);
    elseif strcmp(model.predim, 'exact')
    model = generate_ftt_correctionfull(model);
    end
elseif strcmp(model.pre, 'NL')
    model = generate_ftt_nl(model);
elseif strcmp(model.pre, 'SSM')
    model = generate_ftt_dirt(model);
end

%%
root = pwd;
save([root '/data/' file],'model');
rmpath([rootref '/' model.name])
