function sams = samx0(thetas, N, model)
% given parameters, sample x0
sams = mvnrnd(zeros(model.dimension,1), eye(model.dimension), N)';
end