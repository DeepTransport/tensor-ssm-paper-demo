function u = matvec_prior_Lt(prior, v)
%MATVEC_PRIOR_LT
%
% L'*v
%
% Tiangang Cui, 17/Jan/2014

switch prior.type
    case {'Field'}
        u = prior_cov_lt(prior, v);
    case {'Basis'}
        u = prior.basis'*v;
end

end