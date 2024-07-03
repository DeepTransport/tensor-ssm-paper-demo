function model = setup(model) % prescribed paras

rng(0)
model.m = 1;
model.n = 1;
truepara = model.pre.para;
model.pre.true = truepara;
model.pre.para = NaN(1, 8);
[x1,x2,x3,x4,x5,x6,x7,x8] = true2ftt(truepara', model);
model.theta = [x1;x2;x3;x4;x5;x6;x7;x8];

% ind = [1,2,3];
% model.pre.ind = ind;
% model.d = length(ind);

model.theta = model.theta(model.pre.ind);
model.pre.para = truepara;
model.pre.para(model.pre.ind) = NaN;

if nnz(isnan(model.pre.para)) ~= model.d
    error("model.d should equal to #the unknown parameters in model.pre.para.")
end
rng('shuffle')

end