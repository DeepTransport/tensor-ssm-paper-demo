function [V,d] = prior_cov_eig(prior)
%COV_EIG
%
% Tiangang Cui, 17/Jan/2014

switch prior.cov_type
    case{'GP'}
        [V,D] = eig(prior.C);
    case{'MRF'}
        opts.issym  = 1;
        opts.isreal = 1;
        [V,D] = eigs(@(dv) prior_cov_c(prior, dv), prior.dof, prior.dof-1, 'LA', opts);
end

[d,ind] = sort( diag(D), 'descend' );
V = V(:, ind);

end
