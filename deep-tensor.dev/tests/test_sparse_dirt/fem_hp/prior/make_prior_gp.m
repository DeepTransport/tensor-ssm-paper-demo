function prior = make_prior_gp(mesh, scale, power, sigma)
%MAKE_PRIOR_GP
% makes the Gaussian prior distribution for a given mesh and a correlation length
%
%%%%%%%%%%%%%%%%%%%% input: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mesh:   
% scale:        the tensor for scaling correlation length
% power:        the power term in the kernel
% sigma:        the standard deviation
%
% Tiangang Cui, 12/May/2012

n = size(mesh.nodes,1);
dist2 = zeros(n);

for i   = 1:n
    x = mesh.nodes(i,:);
    d = (mesh.nodes - repmat(x,n,1));
    dist2(:,i)  = sum((d*scale).*d, 2);
end

if power < 0
    nu = abs(power);
    % dist2 = ((pmx-pmx').^2 + (pmy-pmy').^2)/(corr_length^2);
    % scale = 1/(corr_length^2)
    C = (2^(1-nu)/gamma(nu)).*((2*nu*dist2).^(nu/2)).*besselk(nu, sqrt(2*nu*dist2));
    C(dist2==0) = 1;
elseif power == 2
    C = (exp(-0.5*dist2) + 1e-10*eye(n));
else
    C = exp(-0.5*dist2.^(0.5*power));
end

prior.C = C*sigma^2;
prior.L = chol(prior.C, 'lower');
prior.cov_type = 'GP';
prior.type = 'Field';
prior.dof  = n;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%