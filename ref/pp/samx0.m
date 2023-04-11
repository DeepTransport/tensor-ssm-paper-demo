function sams = samx0(thetas, N, model)
% given parameters, sample x0
sams = mvnrnd([50,5]', [25,0;0,1], N)';
end