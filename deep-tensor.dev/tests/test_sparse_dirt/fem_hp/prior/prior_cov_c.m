function u = prior_cov_c(prior, v)
%COV_CV
%
% Tiangang Cui, 17/Jan/2014

switch prior.cov_type
    case {'MRF'}
        t = zeros(size(v));
        u = zeros(size(v));
        t(prior.pk,:) = prior.Lk'\(prior.Lk\v(prior.pk,:));
        t = prior.M'*t;
        u(prior.pk,:) = prior.Lk'\(prior.Lk\t(prior.pk,:));
    case {'GP'}
        u = prior.C*v;
end

end