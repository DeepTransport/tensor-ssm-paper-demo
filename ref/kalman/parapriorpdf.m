function pdf = parapriorpdf(para, model)
% used in SMC^2
% the prior pdf for only parameters
n = size(para,2);
pdf = zeros(1,n);
for k = 1:n
    pdf(k) =  mvnpdf(para(1:model.td,k)');
end
end