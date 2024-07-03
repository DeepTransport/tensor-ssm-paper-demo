function model = setup(model) % prescribed paras

rng(0)
model.d = 2;
model.m = 1;
model.n = 1;
model.theta = [0.6; 0.4];
model.type = 0;

model.theta(1) = norminv((model.theta(1) - 0.1)/0.8);
model.theta(2) = norminv((model.theta(2) - 0.1)/0.8);

rng('shuffle')

end