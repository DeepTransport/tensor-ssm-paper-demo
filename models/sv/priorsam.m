function sams = priorsam(model, N)
sams = [randn(model.d, N) ; mvnrnd(zeros(1,1), 1.25, N)'];
end