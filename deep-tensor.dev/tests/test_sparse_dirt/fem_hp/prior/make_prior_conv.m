function prior = make_prior_conv(mesh, centers, radii, weights, kernel_func)

n  = size(centers, 1);
np = size(mesh.nodes, 1);
prior.basis = zeros(np, n);

for i = 1:n
    x = centers(i,:);
    d = mesh.nodes - repmat(x,np,1);
    s = sqrt(sum(d.^2, 2))/radii(i);
    prior.basis(:,i) = kernel_func(s(:))*weights(i);
end

[Q,R] = qr(prior.basis, 0);
prior.basis_w = (R\Q')';
%
prior.type = 'Basis';
prior.dof  = n;
prior.cov_type = 'Conv';

end

