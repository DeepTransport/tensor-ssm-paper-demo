function model = setup(model) % prescribed paras

rng(0)
model.theta = norminv(([.8; .5]-0.4)/0.6);
model.type = 2;

model.pre.errcov = eye(model.m);
model.pre.C = rand(model.n, model.m);
% model.pre.B = eye(model.m);
% model.pre.C = eye(model.n);%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle')

end